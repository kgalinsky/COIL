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

sub barcodes2numeric {
    my ( $major_alleles, $barcodes ) =
      validate_pos( @_, { type => Params::Validate::ARRAYREF }, $VAL_BARCODES );
    my $converter = _make_converter_obj( $_[0] );
    [ map { _barcode2numeric( $converter, $_ ) } @{ $_[1] } ];
}

sub _make_converter_obj {
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
