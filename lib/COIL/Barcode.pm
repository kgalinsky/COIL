package COIL::Barcode;

use strict;
use warnings;

use Carp;
use List::MoreUtils 'pairwise';
use Params::Validate;

use COIL;
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
        my ($barcode_str) = (split)[-1];
        push @barcodes, [ split //, $barcode_str ];
    }

    croak "Invalid barcodes found"
      unless ( val_barcodes( \@barcodes ) );

    return \@barcodes;
}

=head2 barcodes2numerics

    my $numerics = barcodes2numerics( \@major_alleles, \@barcodes );

Convert barcodes to numeric representations

=cut

sub barcodes2numerics {
    my ( $barcodes, $major_alleles ) =
      validate_pos( @_, $VAL_BARCODES, $VAL_STRAIN );

    unless ( @$major_alleles == @{ $barcodes->[0] } ) {
        croak q{Major allele array doesn't correspond to barcode length.};
    }

    my $struct = _barcode2numeric_struct( $major_alleles );
    [ map { _barcode2numeric( $struct, $_ ) } @$barcodes ];
}

sub _barcode2numeric_struct {
    [
        map {
            my @a;
            @a[ ( 65, 67, 71, 84, 78, 88 ) ] = ( (1) x 4, 2, 3 );
            $a[ ord($_) ] = 0;
            \@a;
        } @{ $_[0] }
    ];
}

sub _barcode2numeric {
    [ pairwise { $a->[ ord($b) ] } @{ $_[0] }, @{ $_[1] } ];
}

1;
