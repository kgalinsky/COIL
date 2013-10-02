#!perl -T
use 5.006;
use strict;
use warnings FATAL => 'all';
use Test::More;

plan tests => 10;

BEGIN {
    use_ok('COIL')           || print "COIL failed\n";
    use_ok('COIL::Validate') || print "COIL::Validate failed\n";
    use_ok('COIL::Barcode')  || print "COIL::Barcode failed\n";
    use_ok('COIL::Likelihood')
      || print "COIL::Likelihood failed\n";
    use_ok('COIL::Pair')        || print "COIL::Pair failed\n";
    use_ok('COIL::Probability') || print "COIL::Probability failed\n";

    use_ok('COIL::Tally::Allele') || print "COIL::Tally::Allele failed\n";
    use_ok('COIL::Likelihood::Allele')
      || print "COIL::Likelihood::Allele failed\n";

    use_ok('COIL::Tally::Pair')      || print "COIL::Tally::Pair failed\n";
    use_ok('COIL::Likelihood::Pair') || print "COIL::Likelihood::Pair failed\n";
}

diag("Testing COIL $COIL::VERSION, Perl $], $^X");
