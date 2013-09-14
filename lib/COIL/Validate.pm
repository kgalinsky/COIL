package COIL::Validate;

use strict;
use warnings;
no warnings 'once';

use Params::Validate;
use Scalar::Util qw/ looks_like_number /;

use Exporter 'import';

our @EXPORT_VAL = qw/
  $VAL_POS_INT
  $VAL_NON_NEG_INT
  $VAL_PROB
  $VAL_POS_REAL

  $VAL_NUC
  $VAL_ALLELE
  $VAL_NALLELE
  $VAL_SALLELE

  $VAL_BARCODE_STR
  $VAL_STRAIN_STR
  $VAL_NUMERIC_STR

  $VAL_BARCODE
  $VAL_STRAIN
  $VAL_NUMERIC

  $VAL_BARCODES
  $VAL_STRAINS
  $VAL_NUMERICS
  /;

our @EXPORT_FUN = qw/
  val_allele
  val_sallele
  val_nallele

  val_barcode
  val_strain
  val_numeric

  val_barcodes
  val_strains
  val_numerics
  /;

our @EXPORT_OK = ( @EXPORT_VAL, @EXPORT_FUN );

our %EXPORT_TAGS = (
    all => \@EXPORT_OK,
    val => \@EXPORT_VAL,
    fun => \@EXPORT_FUN
);

# validate numbers
our $VAL_POS_INT     = { regex => qr/^[1-9]\d*$/ };
our $VAL_NON_NEG_INT = { regex => qr/^\d+$/ };
our $VAL_PROB        = {
    callbacks => {
        'looks like a number' => sub { looks_like_number( $_[0] ) },
        'in [0,1]'            => sub { ( $_[0] >= 0 ) && ( $_[0] <= 1 ) }
    }
};
our $VAL_POS_REAL = {
    callbacks => {
        'looks like a number' => sub { looks_like_number( $_[0] ) },
        '>0'                  => sub { ( $_[0] > 0 ) }
    }
};

# validate strings
our $RE_NUC     = qr/[ACGT]/;
our $RE_SALLELE = qr/^$RE_NUC$/;
our $RE_STRAIN  = qr/^$RE_NUC+$/;

our $RE_UNIT    = qr/[ACGTNX]/;
our $RE_ALLELE  = qr/^$RE_UNIT$/;
our $RE_BARCODE = qr/^$RE_UNIT+$/;

our $RE_NUM     = qr/[0-3]/;
our $RE_NALLELE = qr/^$RE_NUM$/;
our $RE_NUMERIC = qr/^$RE_NUM+$/;

our $VAL_ALLELE  = { regex => $RE_ALLELE };
our $VAL_SALLELE = { regex => $RE_SALLELE };
our $VAL_NALLELE = { regex => $RE_NALLELE };

our $VAL_BARCODE_STR = { regex => $RE_BARCODE };
our $VAL_STRAIN_STR  = { regex => $RE_STRAIN };
our $VAL_NUMERIC_STR = { regex => $RE_NUMERIC };

# validate barcode-type object
sub _gen_VAL_BARCODE {
    my ($val_fun) = @_;

    return {
        type      => Params::Validate::ARRAYREF,
        callbacks => { 'contains valid alleles' => $val_fun }
    };
}

our $VAL_BARCODE = _gen_VAL_BARCODE( \&_val_alleles );
our $VAL_STRAIN  = _gen_VAL_BARCODE( \&_val_salleles );
our $VAL_NUMERIC = _gen_VAL_BARCODE( \&_val_nalleles );

# validate arrays of barcode-type objects
sub _gen_VAL_BARCODES {
    my ( $name, $val_fun ) = @_;
    my $callback_name = "contains valid $name";

    return {
        type      => Params::Validate::ARRAYREF,
        callbacks => {
            $callback_name         => $val_fun,
            'barcodes same length' => \&_val_barcodes_same_length
        }
    };
}

our $VAL_BARCODES = _gen_VAL_BARCODES( 'barcodes', \&_val_barcodes );
our $VAL_STRAINS  = _gen_VAL_BARCODES( 'strains',  \&_val_strains );
our $VAL_NUMERICS = _gen_VAL_BARCODES( 'numerics', \&_val_numerics );

# functions for validating individual alleles
sub val_allele  { $_[0] =~ $RE_ALLELE }
sub val_sallele { $_[0] =~ $RE_SALLELE }
sub val_nallele { $_[0] =~ $RE_NALLELE }

# functions for validating objects
sub _gen_val_barcode {
    my ($val_fun) = @_;

    return sub {
        return unless ( ref( $_[0] ) eq 'ARRAY' );
        return unless &$val_fun( $_[0] );
        return 1;
    };
}

*val_barcode = _gen_val_barcode( \&_val_alleles );
*val_strain  = _gen_val_barcode( \&_val_salleles );
*val_numeric = _gen_val_barcode( \&_val_nalleles );

# functions to validate arrays of barcode-type objects
sub _gen_val_barcodes {
    my ($val_fun) = @_;

    return sub {
        return unless ( ref( $_[0] ) eq 'ARRAY' );
        return unless &$val_fun( $_[0] );
        return unless _val_barcodes_same_length( $_[0] );
        return 1;
    };
}

sub _val_barcodes_same_length {
    my $l = scalar( @{ $_[0][0] } );
    for ( my $i = 1 ; $i < @{ $_[0] } ; $i++ ) {
        return if ( $l != scalar( @{ $_[0][1] } ) );
    }
    return 1;
}

*val_barcodes = _gen_val_barcodes( \&_val_barcodes );
*val_strains  = _gen_val_barcodes( \&_val_strains );
*val_numerics = _gen_val_barcodes( \&_val_numerics );

# hidden functions for validatining individual elements within an array
sub _gen_val_mult {
    my ($val_fun) = @_;

    return sub {
        foreach my $e ( @{ $_[0] } ) { return unless &$val_fun($e) }
        return 1;
    };
}

*_val_alleles  = _gen_val_mult( \&val_allele );
*_val_salleles = _gen_val_mult( \&val_sallele );
*_val_nalleles = _gen_val_mult( \&val_nallele );
*_val_barcodes = _gen_val_mult( \&val_barcode );
*_val_strains  = _gen_val_mult( \&val_strain );
*_val_numerics = _gen_val_mult( \&val_numeric );

1;
