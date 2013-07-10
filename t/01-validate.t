#!perl -T

use strict;
use warnings FATAL => 'all';

use Test::More tests => 17;

use COIL::Validate;

ok( COIL::Validate::val_allele('A') );
ok( COIL::Validate::val_allele('C') );
ok( COIL::Validate::val_allele('G') );
ok( COIL::Validate::val_allele('T') );
ok( COIL::Validate::val_allele('N') );
ok( COIL::Validate::val_allele('X') );
ok( !COIL::Validate::val_allele('Z') );

ok( COIL::Validate::val_barcode( [qw/ A C G T N X /] ) );
ok( !COIL::Validate::val_barcode('scalar') );
ok( !COIL::Validate::val_barcode( [qw/ A C G T N X Z /] ) );

ok( COIL::Validate::val_barcodes( [ [qw/ A C G T N X /] ] ) );
ok(
    COIL::Validate::val_barcodes(
        [ [qw/ A C G T N X /], [qw/ A C G T N X /] ]
    )
);
ok( !COIL::Validate::val_barcodes('scalar') );
ok( !COIL::Validate::val_barcodes( ['scalar'] ) );
ok( !COIL::Validate::val_barcodes( [ [qw/ A C G T N Z /] ] ) );
ok(
    !COIL::Validate::val_barcodes(
        [ [qw/ A C G T N X /], [qw/ A C G T N Z /] ]
    )
);
ok(
    !COIL::Validate::val_barcodes(
        [ [qw/ A C G T N X /], [qw/ A C G T N X X /] ]
    )
);
