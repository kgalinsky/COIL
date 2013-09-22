package COIL::Likelihood;

use strict;
use warnings;

use Params::Validate;
use COIL::Validate ':val';

=head1 NAME

    COIL::Likelihood - base class for likelihood objects

=head1 DESCRIPTION

Likelihood-calculating objects take a barcode in numeric form and calculate
the log likelihood. There are 3 levels of objects in one likelihood object:

    COIL::Likelihood        - calculate l(G|C) at all C
    COIL::Likelihood::Level - calculate l(G|C) at one C
    COIL::Likelihood::Unit  - calculate l(G_i|C)

=cut

=head1 CONSTRUCTOR

=head2 new_from_tally

    my $CLX = COIL::Likelihod::CLASS->new_from_tally( $CTX );
    my $CLX = COIL::Likelihod::CLASS->new_from_tally( $CTX,
                                                      {
                                                          max_COI => 5,
                                                          padding => 0.5
                                                      }
    );


Constructor for COIL::Likelihood::Allele and COIL::Likelihood::Pair. Create a
likelihood calculator from a tally. The following steps are performed:

    Create base level from tally by calling _new_from_tally in the level
    Increment levels until max COI

The padding indicates how much to pad the counts in the tally.

=cut

sub new_from_tally {
    my $class = shift;
    my ( $CTX, @p ) = validate_pos( @_, 1, { default => {} } );
    my %p = validate(
        @p,
        {
            max_COI => { default => 5,   %$VAL_POS_INT },
            padding => { default => 0.5, %$VAL_NON_NEG_REAL }
        }
    );

    my $max_COI = $p{max_COI};
    my $padding = $p{padding};

    my $self = bless [], $class;
    my $L1 = my $L = "${class}::Level"->_new_from_tally( $CTX, $p{padding} );
    push @$self, $L;

    for ( my $i = 1 ; $i < $max_COI ; $i++ ) {
        $L = $L->_increment($L1);
        push @$self, $L;
    }

    return $self;
}

=head1 METHODS

=cut

sub add_error {
    my $self = shift;
    my ($e) = validate_pos(
        @_,
        {
            default => 0.05,
            %$VAL_PROB
        }
    );

    # $E->[$i][$j] = P(G*=i|G=j)
    my $e1 = 1 - $e;
    my $e2 = $e / 2;
    my @E  = ( [ $e1, $e2, $e2 ], [ $e2, $e1, $e2 ], [ $e2, $e2, $e1 ] );

    return bless [ map { $_->_add_error( \@E ) } @$self ], ref($self);
}

=head2 numeric_likelihood

    my $likelihood  = $CL->numeric_likelihood( $numeric );

=head2 numerics_likelihoods

    my $likelihoods = $CL->numerics_likelihoods( \@numerics );

Compute log likelihood for a numeric or an array of numerics.

=cut

sub numeric_likelihood {
    shift->_numeric_likelihood( validate_pos( @_, $VAL_NUMERIC ) );
}

sub numerics_likelihoods {
    my $self = shift;
    my ($numerics) = validate_pos( @_, $VAL_NUMERICS );

    return ( [ map { $self->_numeric_likelihood($_) } @$numerics ] );
}

sub _numeric_likelihood {
    [ map { $_->_numeric_likelihood( $_[1] ) } @{ $_[0] } ];
}

package COIL::Likelihood::Level;

use strict;
use warnings;

use Carp;

use List::MoreUtils 'pairwise';

=head1 NAME

    COIL::Likelihod::Level

=head1 DESCRIPTION

=cut

=head1 CONSTRUCTORS

=cut

=head2 new_from_tally

=cut

sub _new_from_tally {
    my $class = shift;
    ( my $class_CLX = $class ) =~ s/::[^:]+$//;
    my ( $CTX, $padding ) = @_;

    bless [ map { "${class_CLX}::Unit"->_new_from_tally( $_, $padding ) }
          @$CTX ], $class;
}

=head2 increment

=cut

sub _increment {
    my ( $self, $L1 ) = @_;
    no warnings 'once';
    bless [ pairwise { $a->_increment($b) } @$self, @$L1 ], ref($self);
}

=head2 add_error

=cut

sub _add_error {
    my ( $self, $E ) = @_;
    bless [ map { $_->_add_error($E) } @$self ], ref($self);
}

=head1 METHODS

=cut

=head2 numeric_likelihood

=cut

sub _numeric_likelihood {
    croak ref( $_[0] ) . ' does not define _numeric_likelihood';
}

package COIL::Likelihood::Unit;

use strict;
use warnings;

use Carp;

sub _new_from_tally {
    croak "$_[0] does not define _new_from_tally";
}

sub _increment {
    croak ref( $_[0] ) . ' does not define _increment';
}

sub _add_error {
    croak ref( $_[0] ) . '$_[0] does not define _add_error';
}

1;
