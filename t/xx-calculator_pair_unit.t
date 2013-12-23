#!perl -T
use 5.006;
use strict;
use warnings FATAL => 'all';

use Test::More;
use Test::Deep;

plan tests => 4;

use COIL::Calculator::Pair;

my $class = 'COIL::Calculator::Pair::Unit';
my $u1 = $class->_new_from_P( [ [ .64, .16 ], [ .16, .04 ] ] );

cmp_deeply(
    $u1,
    bless(
        [ num( log(.2), .0001 ), num( log(.8), .0001 ), '-inf', 0 ], $class
    ),
    '_new_from_p'
);

my $u2 = $u1->_increment($u1);

cmp_deeply(
    $u2,
    bless(
        [
            num( log(.04), .0001 ),
            num( log(.64), .0001 ),
            num( log(.32), .0001 ),
            0
        ],
        $class
    ),
    '_increment'
);

my $u1e = $u1->_add_error(
    [ [ .95, .01, .04 ], [ .01, .95, .04 ], [ .05, .05, .90 ] ] );

cmp_deeply(
    $u1e,
    bless(
        [
            num( log(.198), .0001 ),
            num( log(.762), .0001 ),
            num( log(.04),  .0001 ),
            0
        ],
        $class
    ),
    '_add_error 1'
);

my $u2e = $u2->_add_error(
    [ [ .95, .01, .04 ], [ .01, .95, .04 ], [ .05, .05, .90 ] ] );

cmp_deeply(
    $u2e,
    bless(
        [
            num( log(.0604), .0001 ),
            num( log(.6244), .0001 ),
            num( log(.3152), .0001 ),
            0
        ],
        $class
    ),
    '_add_error 2'
);

