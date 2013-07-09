package COIL::Barcode;

use strict;
use warnings;

use Carp;

use COIL;

=head1 FUNCTIONS

=head2 validate_barcode

    validate_barcode($barcode);

Validate the barcode and croak if bad.

=head2 validate_barcodes

    validate_barcodes($barcodes);

Validate the barcodes and croak if bad.

=cut

sub validate_barcode {
    croak 'Barcode not arrayref' unless ( ref( $_[0] ) eq 'ARRAY' );

    foreach my $allele ( @{ $_[0] } ) {
        croak qq{Invalid allele "$allele"} unless ( $allele =~ m/^[ACGTNX]$/ );
    }
}

sub validate_barcodes {
    foreach my $barcode ( @{ $_[0] } ) {
        validate_barcode($barcode);
    }
}

=head2 read_barcodes

    my $barcodes = read_barcodes($filename);
    my $barcodes = read_barcodes($filehandle);

Read barcodes from file.

=cut

sub read_barcodes {
    my $fh = COIL::_fh( \@_ );

    my @barcodes;
    while (local $_ = <$fh>) {
        my ($barcode_str) = (split)[-1];
        push @barcodes, [ split //, $barcode_str ];
    }

    validate_barcodes( \@barcodes );
    return \@barcodes;
}

1;
