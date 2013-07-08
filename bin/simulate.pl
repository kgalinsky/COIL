#!/usr/bin/env perl

use strict;
use warnings;

my $ARGV_length = 24;  # barcode length
my $ARGV_MAF    = 0.2; # minor allele frequency
my $ARGV_number = 100; # number of barcodes

for (my $i = 0; $i < $ARGV_number; $i++) {
    for (my $j = 0; $j < $ARGV_length; $j++) {
        print rand() < $ARGV_MAF ? 'A' : 'T';
    }
    print "\n";
}

__END__

=head1 NAME

COIL/simulate.pl - simulate barcodes

=cut
