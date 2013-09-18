#!/usr/bin/env perl

use strict;
use warnings;

use COIL::Barcode;
use COIL::Tally::Allele;
use COIL::Likelihood::Allele;
use COIL::Probability;

use List::MoreUtils 'pairwise';

my $barcodes    = COIL::Barcode::read_barcodes( $ARGV[0] );
my $TA          = COIL::Tally::Allele->tally_barcodes($barcodes);
my $LA          = COIL::Likelihood::Allele->tally2likelihood($TA);
my $LAE         = $LA->add_error();
my $numerics    = $TA->barcodes2numerics($barcodes);
my $likelihoods = $LAE->numerics_likelihoods($numerics);
my $CP          = COIL::Probability->uniform(5);
my @posteriors  = map { $CP->posterior($_) } @$likelihoods;
my @MAPs        = map { $_->mode } @posteriors;

no warnings 'once';
my @Cs = pairwise { $a->credible_interval( 0.95, $b ) } @posteriors, @MAPs;

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
      $MAPs[$i], exp( $posteriors[$i][ $MAPs[$i] - 1 ] ),
      @{ $Cs[$i] };
}
