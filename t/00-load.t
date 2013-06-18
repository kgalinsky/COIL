#!perl -T
use 5.006;
use strict;
use warnings FATAL => 'all';
use Test::More;

plan tests => 2;

BEGIN {
    use_ok('COIL')         || print "COIL failed\n";
    use_ok('COIL::Allele') || print "COIL::Allele failed\n";
}

diag("Testing COIL $COIL::VERSION, Perl $], $^X");
