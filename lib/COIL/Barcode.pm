package COIL::Barcode;

use strict;
use warnings;

use Carp;
use List::MoreUtils 'pairwise';
use Params::Validate;

use COIL '_fh';
use COIL::Validate ':all';

=head1 FUNCTIONS

=head2 read_barcodes

    my $barcodes = read_barcodes($filename);
    my $barcodes = read_barcodes($filehandle);

Read barcodes from file.

=cut

sub read_barcodes {
    my $fh = COIL::_fh( \@_ );

    my @barcodes;
    while ( local $_ = <$fh> ) {
        next if (m/^#/);
        my ($barcode_str) = (split)[-1];
        push @barcodes, [ split //, $barcode_str ];
    }

    croak "Invalid barcodes found"
      unless ( val_barcodes( \@barcodes ) );

    return \@barcodes;
}

1;
