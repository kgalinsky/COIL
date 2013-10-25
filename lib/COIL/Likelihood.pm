package COIL::Likelihood;

use strict;
use warnings;

use Carp;

use Params::Validate;
use COIL::Validate ':val';

=head1 NAME

COIL::Likelihood - base class for likelihood objects

=head1 SYNOPSIS

    my $CLX         = COIL::Likelihood::CLASS->new_from_tally( $CTX );
    my $CLX_E       = $CLX->add_error( $e );
    my $likelihood  = $CLX->numeric_likelihood( $numeric );
    my $likelihoods = $CLX->numerics_likelihoods( $numerics );

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
    my $CLX = COIL::Likelihod::CLASS->new_from_tally(
        $CTX,
        {
            max_COI => 5,
            padding => 1
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
            padding => $VAL_PADDING
        }
    );

    my $max_COI = $p{max_COI};
    my $padding = $p{padding};

    my $self = bless [], $class;
    my $L1 = my $L = "${class}::Level"->_new_from_tally( $CTX, $padding );
    push @$self, $L;

    for ( my $i = 1 ; $i < $max_COI ; $i++ ) {
        $L = $L->_increment($L1);
        push @$self, $L;
    }

    return $self;
}

=head2 add_error

    my $CLX_E = $CLX->add_error();
    my $CLX_E = $CLX->add_error( $e ); # default e=.05
    my $CLX_E = $CLX->add_error( [ $e0, $e+, $e- ]);

Add errors to the likeliood-computing object. The new object is computed by
doing a matrix multiplication (done at the unit level). The matrix may be
specified simply by using one error rate (e) or by specifying 3 error rates.

    e0: ref<->alt conversion
    e+: ref/alt->het conversion
    e-: het->ref/alt conversion (e-)

Let G_i* be the observed allele and A*/a*/N* be the associated events from
before. The matrix is just a stochastic matrix.

    E[i][j] = P(G_i*=j|G_i=i)

        / P(A*|A), P(a*|A), P(N*|A) \
    E = | P(A*|a), P(a*|a), P(N*|a) |
        \ P(A*|N), P(a*|N), P(N*|N) /

        / 1-e, e/2, e/2 \
      = | e/2, 1-e, e/2 |, simple error rate case
        \ e/2, e/2, 1-e /

        / 1-e0-e+, e0,      e+   \
      = | e0     , 1-e0-e+, e+   |, a bit more complicated
        \ e-/2   , e-/2,    1-e- /


=cut

sub add_error {
    shift->_add_error( _val_error( \@_ ) );
}

sub _add_error {
    return bless [ map { $_->_add_error( $_[1] ) } @{ $_[0] } ], ref( $_[0] );

}

# validate errors

sub _val_error {
    my ($p) = @_;
    my $e = shift(@$p);

    unless ( defined($e) ) { $e = 0.05 }

    my $type = ref($e);
    if ($type) {
        if ( $type eq 'ARRAY' ) {
            croak 'Error array must specify e0/e+/e-' unless ( @$e == 3 );
            my ( $e0, $ep, $em ) = @$e;
            return [
                [ 1 - $e0 - $ep, $e0,           $ep ],
                [ $e0,           1 - $e0 - $ep, $ep ],
                [ $em / 2,       $em / 2,       1 - $em ]
            ];
        }
        else { croak qq{"${type}REF" not supported } }
    }
    else {
        my $e1 = 1 - $e;
        my $e2 = $e / 2;
        return [ [ $e1, $e2, $e2 ], [ $e2, $e1, $e2 ], [ $e2, $e2, $e1 ] ];
    }
}

=head1 METHODS

=cut

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

=head1 NAME

COIL::Likelihood::Level - base class for likelihood levels

=cut

package COIL::Likelihood::Level;

use strict;
use warnings;

use Carp;

use List::MoreUtils 'pairwise';

=head1 DESCRIPTION

This module does the actual likelihood computation for a numeric. It is
composed of multiple units. Most of the constructors/methods in this module
actually organize calls to the Unit class.

=cut

=head1 CONSTRUCTORS

=cut

=head2 new_from_tally

    COIL::Likelihood::CLASS::Level->new_from_tally( $tally );

Create a new level.

=cut

sub new_from_tally {
    my $class = shift;
    ( my $tally_class = $class ) =~
      s/COIL::Likelihood::([^:]+)::Level/COIL::Tally::$1/;
    $class->_new_from_tally( validate_pos( @_, { isa => $tally_class } ) );
}

sub _new_from_tally {
    my $class = shift;
    ( my $unit_class = $class ) =~ s/Level$/Unit/;
    my ( $CTX, $padding ) = @_;

    bless [ map { $unit_class->_new_from_tally( $_, $padding ) } @$CTX ],
      $class;
}

=head2 increment

Increment a level by zipping the current level and the base level.

=cut

sub increment {
    my $self = shift;
    $self->_increment( validate_pos( @_, { isa => ref($self) } ) );
}

sub _increment {
    my ( $self, $L1 ) = @_;
    no warnings 'once';
    bless [ pairwise { $a->_increment($b) } @$self, @$L1 ], ref($self);
}

=head2 add_error

Create a new level by adding errors to the current level.

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

{
    no warnings 'once';
    *add_error          = \&COIL::Likelihood::add_error;
    *numeric_likelihood = \&COIL::Likelihood::numeric_likelihood;
}

=head1 NAME

COIL::Likelihood::Unit - base class for likelihood units

=cut

package COIL::Likelihood::Unit;

use strict;
use warnings;

use Params::Validate;

use Carp;

=head1 DESCRIPTION

Likelihoods for individual genotype calls at specific COIs are stored in a
unit object. Computing a likelihood for a numeric involves doing array lookups.
This module actually handles the intial computations, which are all
constructors. The base class will try handle as much as it can as far as
parameter validations, but the following constructors must be defined in the
sub-classes:

    _new_from_tally
    _increment
    _add_error

=cut

=head1 CONSTRUCTORS

=cut

=head2 new_from_tally

    COIL::Likelihood::CLASS::Unit->new_from_tally( $tally_unit );

Create a COIL::Likelihood::CLASS::Unit from a COIL::Tally::CLASS::Unit object.
The CLASS must match.

=cut

sub new_from_tally {
    my $class = shift;
    ( my $tally_class = $class ) =~ s/Likelihood/Tally/;
    $class->_new_from_tally( validate_pos( @_, { isa => $tally_class } ) );
}

sub _new_from_tally {
    croak "$_[0] does not define _new_from_tally";
}

=head2 increment

    my $U_n   = $U_n_1->increment( $U_1 );
    my $U_a_b = $U_a->increment( $U_b );

Create a unit for COI=a+b from units where COI=a and COI=b. The primary use is
to create a likelihood object iterating through each COI.

=cut

sub _increment {
    croak ref( $_[0] ) . ' does not define _increment';
}

sub _add_error {
    croak ref( $_[0] ) . '$_[0] does not define _add_error';
}

{
    no warnings 'once';
    *increment = \&COIL::Likelihood::Level::increment;
    *add_error = \&COIL::Likelihood::add_error;
}

1;
