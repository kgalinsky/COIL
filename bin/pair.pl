#!/usr/bin/env perl

use strict;
use warnings;

use List::MoreUtils 'pairwise';

use COIL::Barcode;
use COIL::Tally::Allele;
use COIL::Tally::Pair;
use COIL::Likelihood::Pair;
use COIL::Probability;

my $barcodes    = COIL::Barcode::read_barcodes( $ARGV[0] );
my $CTA         = COIL::Tally::Allele->new_from_barcodes($barcodes);
my $numerics    = $CTA->barcodes2numerics($barcodes);
my $CTP         = COIL::Tally::Pair->new_from_numerics($numerics);
my $CLP         = COIL::Likelihood::Pair->new_from_tally($CTP);
my $CLP_E       = $CLP->add_error();
my $likelihoods = $CLP_E->numerics_likelihoods($numerics);
my $CP          = COIL::Probability->uniform(5);
my @posteriors  = map { $CP->posterior($_) } @$likelihoods;
my @MAPs        = map { $_->mode } @posteriors;
no warnings 'once';
my @Cs = pairwise { $a->credible_interval( 0.95, $b ) } @posteriors, @MAPs;

local $, = "\t";
local $\ = "\n";

for ( my $i = 0 ; $i < @$barcodes ; $i++ ) {
    print
      join( '', @{ $barcodes->[$i] } ),
      join( '', @{ $numerics->[$i] } ),
      $MAPs[$i], exp( $posteriors[$i][ $MAPs[$i] - 1 ] ),
      @{ $Cs[$i] };
}
