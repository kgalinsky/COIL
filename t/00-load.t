#!perl -T
use 5.006;
use strict;
use warnings FATAL => 'all';
use Test::More;

plan tests => 10;

sub use_or_fail {
    use_ok( $_[0] ) || print "$_[0] failed\n";
}

BEGIN {
    foreach my $p (
        qw/
        COIL

        COIL::Validate
        COIL::Barcode
        COIL::Calculator
        COIL::Pair
        COIL::Probability

        COIL::Tally::Allele
        COIL::Calculator::Allele

        COIL::Tally::Pair
        COIL::Calculator::Pair
        /
      )
    {
        use_or_fail($p);
    }
}

diag("Testing COIL $COIL::VERSION, Perl $], $^X");
