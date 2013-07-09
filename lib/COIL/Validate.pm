package COIL::Validate;

use strict;
use warnings;

use Exporter 'import';
our @EXPORT_OK = qw/
  $RE_UNIT
  $RE_ALLELE
  $RE_BARCODE

  $VAL_ALLELE
  $VAL_BARCODE_STR
  $VAL_BARCODE
  $VAL_BARCODES

  val_allele
  val_barcode
  val_barcodes
  /;

our %EXPORT_TAGS = (
    all => \@EXPORT_OK,
    re  => [
        qw/
          $RE_UNIT
          $RE_ALLELE
          $RE_BARCODE
          /
    ],
    val => [
        qw/
          $VAL_ALLELE
          $VAL_BARCODE_STR
          $VAL_BARCODE
          $VAL_BARCODES
          /
    ],
    fun => [
        qw/
          val_allele
          val_barcode
          val_barcodes
          /
    ]
);

use Params::Validate;

our $RE_UNIT    = qr/[ACGTNX]/;
our $RE_ALLELE  = qr/^$RE_UNIT$/;
our $RE_BARCODE = qr/^$RE_UNIT+$/;

our $VAL_ALLELE = { type => Params::Validate::SCALAR, regex => $RE_ALLELE };
our $VAL_BARCODE_STR =
  { type => Params::Validate::SCALAR, regex => $RE_BARCODE };

our $VAL_BARCODE = {
    type      => Params::Validate::ARRAYREF,
    callbacks => {
        'allele validation' =>
          sub { validate_pos( @{ $_[0] }, ($VAL_ALLELE) x @{ $_[0] } ) }
    }
};

our $VAL_BARCODES = {
    type      => Params::Validate::ARRAYREF,
    callbacks => {
        'barcode validation' =>
          sub { validate_pos( @{ $_[0] }, ($VAL_BARCODE) x @{ $_[0] } ) },
        'barcodes same length' => sub {
            my %l;
            foreach my $b ( @{ $_[0] } ) { $l{ scalar(@$b) } = 1 }
            keys(%l) == 1;
          }
    }
};

sub val_allele   { validate_pos( @_, $VAL_ALLELE ) }
sub val_barcode  { validate_pos( @_, $VAL_BARCODE ) }
sub val_barcodes { validate_pos( @_, $VAL_BARCODES ) }

1;
