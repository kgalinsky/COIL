#!perl -T

use strict;
use warnings FATAL => 'all';

use Test::More tests => 34;

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

ok( COIL::Validate::val_sallele('A') );
ok( COIL::Validate::val_sallele('C') );
ok( COIL::Validate::val_sallele('G') );
ok( COIL::Validate::val_sallele('T') );
ok( !COIL::Validate::val_sallele('N') );
ok( !COIL::Validate::val_sallele('X') );
ok( !COIL::Validate::val_sallele('Z') );

ok( COIL::Validate::val_strain( [qw/ A C G T /] ) );
ok( !COIL::Validate::val_strain('scalar') );
ok( !COIL::Validate::val_strain( [qw/ A C G T N X Z /] ) );

ok( COIL::Validate::val_strains( [ [qw/ A C G T /] ] ) );
ok(
    COIL::Validate::val_strains(
        [ [qw/ A C G T /], [qw/ A C G A /] ]
    )
);
ok( !COIL::Validate::val_strains('scalar') );
ok( !COIL::Validate::val_strains( ['scalar'] ) );
ok( !COIL::Validate::val_strains( [ [qw/ A C G T N Z /] ] ) );
ok(
    !COIL::Validate::val_strains(
        [ [qw/ A C G T /], [qw/ A C G Z /] ]
    )
);
ok(
    !COIL::Validate::val_strains(
        [ [qw/ A C G T /], [qw/ A C G T T /] ]
    )
);