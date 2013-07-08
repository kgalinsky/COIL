package COIL::Barcode;

use strict;
use warnings;

use Carp;

=head1 FUNCTIONS

=head2 validate_barcode

    validate_barcode($barcode);

Validate the barcode and croak if bad.

=head2 validate_barcodes

    validate_barcodes($barcodes);

Validate the barcodes and croak if bad.

=cut

sub validate_barcode {
    croak 'Barcodes are arrayrefs' unless ( ref( $_[0] ) eq 'ARRAY' );

    foreach my $allele ( @{ $_[0] } ) {
        croak qq{Invalid allele "$allele"} unless ( $allele =~ m/^[ACGTNX]$/ );
    }
}

sub validate_barcodes {
    foreach my $barcode ( @{ $_[0] } ) {
        validate_barcode($barcode);
    }
}

1;
