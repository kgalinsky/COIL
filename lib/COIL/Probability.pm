package COIL::Probability;

use strict;
use warnings;

use Carp;

use List::Util qw/ max min reduce sum /;
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
            max_COI => 10,    # hard limit
            CDF     => .99,   # CDF >= .99
            PDF     => .001   # PDF <  .001
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
            max_COI => { optional => 1, %$VAL_POS_INT },
            CDF     => { optional => 1, %$VAL_PROB },
            PDF     => { optional => 1, %$VAL_PROB }
        }
    );

    $p{max_COI} = 5 unless ( $p{max} || $p{CDF} || $p{PDF} );
    $p{max_COI} ||= 100;
    $p{CDF}     ||= 1;
    $p{PDF}     ||= 0;

    # p = f(x)/(1-f(0))
    # log(p) = log(f(x)) - log(1-f(0))
    #        = k*log(lambda) - lambda - log(x!) - log(1-f(0))
    my $log_p      = -1 * $lambda - log( 1 - exp( -1 * $lambda ) );
    my $CDF        = 0;
    my $log_lambda = log($lambda);

    my $self = bless [], $class;
    for ( my $i = 1 ; $i <= $p{max_COI} ; $i++ ) {
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

    my $CP = COIL::Probability->uniform( $max_COI );

Discrete Uniform( 1, max_COI ).

=cut

sub uniform {
    my $class = shift;
    my ($max_COI) = validate_pos( @_, $VAL_POS_INT );
    return bless [ ( log( 1 / $max_COI ) ) x $max_COI ], $class;
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

    # compute posterior
    # P(A|B) ~ L(A|B)P(A)
    # log P(A|B) ~ log(L(A|B)) + log(P(A))
    no warnings 'once';
    my $posterior = [ pairwise { $a + $b } @$self, @$l ];

    # first scale so max log prob is 0 to avoid rounding errors
    my $scale = max @$posterior;
    $_ -= $scale foreach (@$posterior);

    # scale so probabilities sum to 1
    $scale = log sum map { exp $_ } @$posterior;
    $_ -= $scale foreach (@$posterior);

    return bless $posterior, ref($self);
}

=head2 mode

    my $COI = $CP->mode();

Return most likely COI.

=cut

sub mode {
    my $self = shift;
    ( reduce { $self->[$a] < $self->[$b] ? $b : $a } ( 0 .. $#$self ) ) + 1;
}

=head2 credible_interval

    my $COI_credible_interval = $CP->credible_interval();
    my $COI_credible_interval = $CP->credible_interval( 0.95 );
    my $COI_credible_interval = $CP->credible_interval( 0.95, $mode );

Return credible interval.

=cut

sub credible_interval {
    my $self = shift;
    my ( $threshold, $mode ) = validate_pos(
        @_,
        { default  => 0.95, %$VAL_PROB },
        { optional => 1,    %$VAL_POS_INT }
    );

    if ($mode) { croak 'Mode out of range' unless ( $mode <= @$self ) }
    else       { $mode = $self->mode() }

    $mode--;
    my $lower = my $upper = $mode;
    my $conf = exp( $self->[$mode] );

    while ( $conf < $threshold ) {
        if    ( $lower == 0 )       { $conf += exp( $self->[ ++$upper ] ) }
        elsif ( $upper == $#$self ) { $conf += exp( $self->[ --$lower ] ) }
        else {
            if ( $self->[ $lower - 1 ] >= $self->[ $upper + 1 ] ) {
                $conf += exp( $self->[ --$lower ] );
            }

            else { $conf += exp( $self->[ ++$upper ] ) }
        }
    }

    [ $lower + 1, $upper + 1, $conf ];
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
