package COIL::Allele;

use strict;
use warnings FATAL => 'all';

use Params::Validate;

sub tally {
    my ($barcodes) = validate_pos( @_, 1 );
    my @tallies;
    foreach my $barcode (@$barcodes) {
        for ( my $i = 0 ; $i < @$barcode ; $i++ ) {
            $tallies[$i]{ $barcode->[$i] }++;
        }
    }
    return \@tallies;
}

1;
