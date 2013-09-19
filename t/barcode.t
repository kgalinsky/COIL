#!perl -T

use strict;
use warnings FATAL => 'all';

use Test::More tests => 1;

use COIL::Barcode;

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

is_deeply( COIL::Barcode::read_barcodes( \$good_barcodes ),
    \@barcodes, 'Read barcodes from scalar ref' );
