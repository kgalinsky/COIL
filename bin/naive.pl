#!/usr/bin/env perl

use strict;
use warnings;

use COIL::Barcode;
use COIL::Allele;
use COIL::Drivers;

my $barcodes = COIL::Barcode::read_barcodes( $ARGV[0] );
my $tallies  = COIL::Allele::tally($barcodes);
my $A        = COIL::Allele->tallies2obj($tallies);

$A->write(*STDERR);

my $lGis     = $A->log_likelihoods(5);
my $lOis     = COIL::Allele::error_log_likelihoods( $lGis, .01 );
my $numerics = COIL::Barcode::barcodes2numerics( $barcodes, $A->majors );
my $lCGs = COIL::Drivers::naive( $lOis, $numerics );    # l(C|G) = log P(G|C)
my $PCGs = COIL::Drivers::posteriors($lCGs);            # P(C|G)
my $MAPs = COIL::Drivers::max_is($PCGs);
my $Cs = COIL::Drivers::credible_intervals( $PCGs, $MAPs );

local $, = "\t";
local $\ = "\n";

for ( my $i = 0 ; $i < @$barcodes ; $i++ ) {
    print
      join( '', @{ $barcodes->[$i] } ),
      join( '', @{ $numerics->[$i] } ),
      $MAPs->[$i] + 1,
      $Cs->[$i][0] + 1,
      $Cs->[$i][1] + 1,
      $Cs->[$i][2];
}
