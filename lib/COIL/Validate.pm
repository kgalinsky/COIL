package COIL::Validate;

use strict;
use warnings;

use Exporter 'import';
our @EXPORT_OK = qw/
  $RE_UNIT
  $RE_ALLELE
  $RE_BARCODE
  $RE_NUC
  $RE_SALLELE

  $VAL_ALLELE
  $VAL_BARCODE_STR
  $VAL_BARCODE
  $VAL_BARCODES
  $VAL_SALLELE
  $VAL_STRAIN
  $VAL_STRAINS

  val_allele
  val_barcode
  val_barcodes
  val_sallele
  val_strain
  val_strains
  /;

our %EXPORT_TAGS = (
    all => \@EXPORT_OK,
    re  => [
        qw/
          $RE_UNIT
          $RE_ALLELE
          $RE_BARCODE
          $RE_NUC
          $RE_SALLELE
          /
    ],
    val => [
        qw/
          $VAL_ALLELE
          $VAL_BARCODE_STR
          $VAL_BARCODE
          $VAL_BARCODES
          $VAL_SALLELE
          $VAL_STRAIN
          $VAL_STRAINS
          /
    ],
    fun => [
        qw/
          val_allele
          val_barcode
          val_barcodes
          val_sallele
          val_strain
          val_strains
          /
    ]
);

use Params::Validate;

our $RE_UNIT    = qr/[ACGTNX]/;
our $RE_ALLELE  = qr/^$RE_UNIT$/;
our $RE_BARCODE = qr/^$RE_UNIT+$/;
our $RE_NUC     = qr/[ACGT]/;
our $RE_SALLELE = qr/^$RE_NUC$/;

our $VAL_ALLELE = { type => Params::Validate::SCALAR, regex => $RE_ALLELE };
our $VAL_BARCODE_STR =
  { type => Params::Validate::SCALAR, regex => $RE_BARCODE };

our $VAL_BARCODE = {
    type      => Params::Validate::ARRAYREF,
    callbacks => { 'contains valid alleles' => \&_val_alleles }
};

our $VAL_BARCODES = {
    type      => Params::Validate::ARRAYREF,
    callbacks => {
        'contains valid barcodes' => \&_val_barcodes,
        'barcodes same length'    => \&_val_barcodes_same_length
    }
};

our $VAL_SALLELE = { type => Params::Validate::SCALAR, regex => $RE_SALLELE };

our $VAL_STRAIN = {
    type      => Params::Validate::ARRAYREF,
    callbacks => { 'contains valid alleles' => \&_val_salleles }
};

our $VAL_STRAINS = {
    type      => Params::Validate::ARRAYREF,
    callbacks => {
        'contains valid barcodes' => \&_val_strains,
        'barcodes same length'    => \&_val_barcodes_same_length
    }
};

sub val_allele { $_[0] =~ $RE_ALLELE }

sub _val_alleles {
    foreach my $a ( @{ $_[0] } ) { return unless val_allele($a) }
    return 1;
}

sub val_barcode {
    return unless ( ref( $_[0] ) eq 'ARRAY' );
    return unless _val_alleles(@_);
    return 1;
}

sub _val_barcodes {
    foreach my $b ( @{ $_[0] } ) { return unless val_barcode($b) }
    return 1;
}

sub val_barcodes {
    return unless ( ref( $_[0] ) eq 'ARRAY' );
    return unless _val_barcodes(@_);
    return unless _val_barcodes_same_length(@_);
    return 1;
}

# compare length of each barcode to that of first
sub _val_barcodes_same_length {
    my $l = scalar( @{ $_[0][0] } );
    for ( my $i = 1 ; $i < @{ $_[0] } ; $i++ ) {
        return if ( $l != scalar( @{ $_[0][1] } ) );
    }
    return 1;
}

sub val_sallele { $_[0] =~ $RE_SALLELE }

sub _val_salleles {
    foreach my $a ( @{ $_[0] } ) { return unless val_sallele($a) }
    return 1;
}

sub val_strain {
    return unless ( ref( $_[0] ) eq 'ARRAY' );
    return unless _val_salleles(@_);
    return 1;
}

sub _val_strains {
    foreach my $b ( @{ $_[0] } ) { return unless val_strain($b) }
    return 1;
}

sub val_strains {
    return unless ( ref( $_[0] ) eq 'ARRAY' );
    return unless _val_strains(@_);
    return unless _val_barcodes_same_length(@_);
    return 1;
}

1;
