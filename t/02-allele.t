#!perl -T

use strict;
use warnings FATAL => 'all';

use Test::More tests => 1;

use COIL::Allele;

my @barcodes = (
    [qw/ A C G /], [qw/ A C G /], [qw/ A C T /], [qw/ A C T /],
    [qw/ A G X /], [qw/ A G X /], [qw/ A G N /], [qw/ A G N /]
);
my @expected_tally =
  ( { A => 8 }, { C => 4, G => 4 }, { G => 2, T => 2, X => 2, N => 2 } );
is_deeply(COIL::Allele::tally(\@barcodes), \@expected_tally, 'tally');
