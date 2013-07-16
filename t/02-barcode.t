#!perl -T

use strict;
use warnings FATAL => 'all';

use Test::More tests => 2;

use COIL::Barcode;

my @major_alleles = qw/ A C G T /;

my $good_barcodes = <<ENDL;
TGGC
ACAC
TNAT
NGGN
NNAC
AGXT
ENDL

my @barcodes = (
    [qw/ T G G C /], [qw/ A C A C /], [qw/ T N A T /], [qw/ N G G N /],
    [qw/ N N A C /], [qw/ A G X T /]
);

my @numerics = (
    [qw/ 1 1 0 1 /], [qw/ 0 0 1 1 /], [qw/ 1 2 1 0 /], [qw/ 2 1 0 2 /],
    [qw/ 2 2 1 1 /], [qw/ 0 1 3 0 /]
);

is_deeply( COIL::Barcode::read_barcodes( \$good_barcodes ),
    \@barcodes, 'Read barcodes from scalar ref' );

is_deeply( COIL::Barcode::barcodes2numerics( \@barcodes, \@major_alleles ),
    \@numerics, 'Convert barcodes to numerics' );