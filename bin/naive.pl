#!/usr/bin/env perl

use strict;
use warnings;

use COIL::Barcode;
use COIL::Tally::Allele;
use COIL::Likelihood::Allele;
use COIL::Likelihood;

my $barcodes = COIL::Barcode::read_barcodes( $ARGV[0] );
my $TA       = COIL::Tally::Allele->tally_barcodes($barcodes);
my $LA       = COIL::Likelihood::Allele->tally2likelihood($TA);
my $LAE        = $LA->add_error();
my $numerics    = $TA->barcodes2numerics($barcodes);
my $likelihoods = $LA->numerics_likelihoods($numerics);
my $posteriors  = COIL::Likelihood::posteriors($likelihoods);
my $MAPs = COIL::Likelihood::MAPs($posteriors);
my $Cs = COIL::Likelihood::credible_intervals( $posteriors, $MAPs );

print STDERR "Tally:\n";
$TA->write(*STDERR);
print STDERR "\n";

print STDERR "Likelihood:\n";
$LA->write(*STDERR);
print STDERR "\n";

print STDERR "Error Likelihood:\n";
$LAE->write(*STDERR);
print STDERR "\n";

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
