package COIL::Calculator;

use strict;
use warnings;

use Carp;

use Params::Validate;
use COIL::Validate ':val';

=head1 NAME

COIL::Calculator - base class for likelihood calculator objects

=head1 SYNOPSIS

    my $calc         = COIL::Calculator::CLASS->new_from_tally( $CTX );
    my $calc_err       = $calc->add_error( $e );
    my $likelihood  = $calc->numeric_likelihood( $numeric );
    my $likelihoods = $calc->numerics_likelihoods( $numerics );

=head1 DESCRIPTION

Likelihood-calculating objects take a barcode in numeric form and calculate
the log likelihood. There are 3 levels of objects in one likelihood object:

    COIL::Calculator        - calculate l(G|C) at all C
    COIL::Calculator::Level - calculate l(G|C) at one C
    COIL::Calculator::Unit  - calculate l(G_i|C) at one C

=cut

=head1 CONSTRUCTOR

=head2 new_from_tally

    my $calc = COIL::Calculator::CLASS->new_from_tally( $CTX );
    my $calc = COIL::Calculator::CLASS->new_from_tally(
        $CTX,
        {
            max_COI => 5,
            padding => 1
        }
    );


Create a likelihood calculator from a tally. Works by creating the base level
for the calculator by creating a calculator unit for each unit in the tally,
then incrementing the levels until the max COI level.

=cut

sub new_from_tally {
    my $class = shift;
    my ( $CTX, @p ) = validate_pos( @_, 1, { default => {} } );
    my %p = validate(
        @p,
        {
            max_COI => { default => 5, %$VAL_POS_INT },
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

    my $calc_err = $calc->add_error();
    my $calc_err = $calc->add_error( $e ); # default e=.05
    my $calc_err = $calc->add_error( [ $e0, $e+, $e- ]);

Add errors to the likelihood-calculating object. The error rate may be
specified as a single error rate (e; default is 5%) or using a 3 error rates:

    e0: ref<->alt conversion
    e+: ref/alt->het conversion
    e-: het->ref/alt conversion (e-)

This method must be implemented in the COIL::Calculator::CLASS::Unit. See
COIL::Calculator::Allele for one documented implementation.

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

    my $likelihood  = $calc->numeric_likelihood( $numeric );

=head2 numerics_likelihoods

    my $likelihoods = $calc->numerics_likelihoods( \@numerics );

Calculate log likelihood for a numeric or an array of numerics.

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

COIL::Calculator::Level - base class for calculator levels

=cut

package COIL::Calculator::Level;

use strict;
use warnings;

use Carp;

use List::MoreUtils 'pairwise';

=head1 DESCRIPTION

This module does the actual likelihood calculatation for a numeric. It is
composed of multiple units. Most of the constructors/methods in this module
actually organize calls to the Unit class/objects.

=cut

=head1 CONSTRUCTORS

=cut

=head2 new_from_tally

    my $level_0 = COIL::Calculator::CLASS::Level->new_from_tally( $tally );

Create the base level from a tally.

=cut

sub new_from_tally {
    my $class = shift;
    ( my $tally_class = $class ) =~
      s/COIL::Calculator::([^:]+)::Level/COIL::Tally::$1/;
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

    my $level_n = $level_n_1->increment( $level_1 );

Increment a level by adding the base level to the current level.

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

    my $level_err = $level->add_error

Create a new level by adding errors to the current level.

=cut

sub _add_error {
    my ( $self, $E ) = @_;
    bless [ map { $_->_add_error($E) } @$self ], ref($self);
}

=head1 METHODS

=cut

=head2 numeric_likelihood

    my $l = $level->numeric_likelihood( $numeric );

Calculate the likelihood of the level's COI for a particular barcode.

=cut

sub _numeric_likelihood {
    croak ref( $_[0] ) . ' does not define _numeric_likelihood';
}

{
    no warnings 'once';
    *add_error          = \&COIL::Calculator::add_error;
    *numeric_likelihood = \&COIL::Calculator::numeric_likelihood;
}

=head1 NAME

COIL::Calculator::Unit - base class for likelihood units

=cut

package COIL::Calculator::Unit;

use strict;
use warnings;

use Params::Validate;

use Carp;

=head1 DESCRIPTION

Likelihoods for individual genotype calls at specific COIs are stored in a
unit object. Calculating a likelihood for a numeric involves doing array
lookups. This module actually handles the intial calculations, which are all
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

    my $unit_1 = COIL::Calculator::CLASS::Unit->new_from_tally( $tally_unit );

Create a COIL::Calculator::CLASS::Unit from a COIL::Tally::CLASS::Unit object.
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

    my $unit_n   = $unit_n_1->increment( $unit_1 );
    my $unit_a_b = $unit_a->increment( $unit_b );

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
    *increment = \&COIL::Calculator::Level::increment;
    *add_error = \&COIL::Calculator::add_error;
}

1;
