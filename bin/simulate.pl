#!/usr/bin/env perl

use strict;
use warnings;

use COIL::Tally::Allele;
use COIL::Likelihood::Allele;

my $ARGV_n1     = 24;      # barcode length
my $ARGV_alpha  = 1;       # alpha for MAF beta distribution
my $ARGV_beta   = 1;       # beta for MAF beta distribution
my $ARGV_n2     = 100;     # number of barcodes
my $ARGV_lambda = 1;       # lambda for COI Poisson distribution
my $ARGV_f      = 0.05;    # failure

my @poisson = ( exp( -1 * $ARGV_lambda ) );
my @n_COI   = (0);

my $c = 1;
while (1) {
    $poisson[$c] = $poisson[ $c - 1 ] * $ARGV_lambda / $c;
    my $n_COI = int( $poisson[$c] / ( 1 - $poisson[0] ) * $ARGV_n2 + .5 );
    last unless ($n_COI);
    push @n_COI, $n_COI;
    $c++;
}

my $CTA =
  COIL::Tally::Allele->random_tally( $ARGV_n1, $ARGV_alpha, $ARGV_beta );
my $CLA =
  COIL::Likelihood::Allele->tally2likelihood( $CTA, { max_COI => $#n_COI } );

local $, = "\t";
local $\ = "\n";

my $i = 1;
for ( my $c = 1 ; $c < @n_COI ; $c++ ) {
    my $numerics =
      [ map { $CLA->random_numeric( $c, $ARGV_f ) } (0) x $n_COI[$c] ];
    my $barcodes = $CTA->numerics2barcodes($numerics);
    print( 'sim' .  $i++, $c, join( '', @$_ ) ) foreach @$barcodes;
}

__END__

=head1 NAME

COIL/simulate.pl - simulate barcodes

=cut
