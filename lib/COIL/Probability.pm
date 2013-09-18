package COIL::Probability;

use strict;
use warnings;

use Carp;

use List::Util qw/ sum /;
use List::MoreUtils qw/ pairwise /;

use Params::Validate;
use COIL::Validate ':val';

=head1 NAME

COIL::Probability - COI probability distributions

=head1 SYNOPSIS

Create and manipulate COI probabilities.

=cut

=head1 CONSTRUCTORS

=cut

=head2 poisson

    my $CP = COIL::Probability->poisson( $lambda );
    my $CP = COIL::Probability->poisson( $lambda,
        { # how we determine the max COI to generate
            max => 10,    # hard limit
            CDF => .99,   # CDF >= .99
            PDF => .001   # PDF <  .001
        } );

Generate a distribution of COIs using a *truncated* Poisson distribution. COI=0
is discarded and COI > $threshold is also discarded. By default, COIs up to 5
are generated, and a hard limit of COI=100 is used (e.g. if CDF=1 or PDF=0 is
used as a threshold).

=cut

sub poisson {
    my $class = shift;
    my ( $lambda, @p ) = validate_pos( @_, $VAL_POS_REAL, { default => {} } );
    my %p = validate(
        @p,
        {
            max => { optional => 1, %$VAL_POS_INT },
            CDF => { optional => 1, %$VAL_PROB },
            PDF => { optional => 1, %$VAL_PROB }
        }
    );

    $p{max} = 5 unless ( $p{max} || $p{CDF} || $p{PDF} );
    $p{max} ||= 100;
    $p{CDF} ||= 1;
    $p{PDF} ||= 0;

    # p = f(x)/(1-f(0))
    # log(p) = log(f(x)) - log(1-f(0))
    #        = k*log(lambda) - lambda - log(x!) - log(1-f(0))
    my $log_p      = -1 * $lambda - log( 1 - exp( -1 * $lambda ) );
    my $CDF        = 0;
    my $log_lambda = log($lambda);

    my $self = bless [], $class;
    for ( my $i = 1 ; $i <= $p{max} ; $i++ ) {
        $log_p += $log_lambda - log($i);
        my $p = exp($log_p);

        last if ( $p < $p{PDF} );

        push @$self, $log_p;
        $CDF += $p;

        last if ( $CDF >= $p{CDF} );
    }

    return $self;
}

=head2 uniform

    my $CP = COIL::Probability->uniform( $max );

Discrete Uniform( 1, max ).

=cut

sub uniform {
    my $class = shift;
    my ($max) = validate_pos( @_, $VAL_POS_INT );
    return bless [ ( log( 1 / $max ) ) x $max ], $class;
}

=head1 METHODS

=cut

=head2 posterior

    my $CP2 = $CP->posterior( $log_likelihood );

=cut

sub posterior {
    my $self = shift;
    my ($l) = validate_pos( @_, $VAL_NON_POS_REALS );

    croak '|posterior| != |likelihood|'
      unless ( $#$self == $#$l );

    no warnings 'once';
    my $posterior = [ pairwise { $a + $b } @$self, @$l ];
    my $scale = log sum map { exp $_ } @$posterior;
    $_ -= $scale foreach (@$posterior);

    return bless $posterior, ref($self);
}

=head2 COIs

    my $COIs = $CP->COIs( $n );

Create a list of almost n COIs where the frequency of each COI is proportional to its
density in the prior. Useful for simulations.

=cut

sub COIs {
    my $self = shift;
    my ($n) = validate_pos( @_, $VAL_POS_INT );
    return [ map { ( $_ + 1 ) x ( exp( $self->[$_] ) * $n ) }
          ( 0 .. $#$self ) ];
}

1;
