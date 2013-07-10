package COIL::Allele;

use strict;
use warnings;

use Params::Validate;
use COIL::Validate ':val';

sub tally {
    my ($barcodes) = validate_pos( @_, $VAL_BARCODES );

    my @tallies;
    foreach my $barcode (@$barcodes) {
        for ( my $i = 0 ; $i < @$barcode ; $i++ ) {
            $tallies[$i]{ $barcode->[$i] }++;
        }
    }

    return \@tallies;
}

1;
