package COIL::Barcode;

use strict;
use warnings;

use Carp;

use Params::Validate;
use COIL::Validate ':all';

use COIL '_fh';

=head1 FUNCTIONS

=head2 read_barcodes

    my $barcodes = read_barcodes($filename);
    my $barcodes = read_barcodes($filehandle);

Read barcodes from file.

=cut

sub read_barcodes {
    my $fh = COIL::_fh( \@_ );

    # kludge to read files with just \r line endings

    # read in entire file
    local $/;
    local $_ = <$fh>;

    my @barcodes =    # read steps backwards
      map  { [ split m// ] }    # split into array of chars
      map  { (split)[-1] }      # grab last column
      grep { !/^#/ }            # remove comments
      split /[\r\n]/;           # split into lines

    croak "Invalid barcodes found"
      unless ( val_barcodes( \@barcodes ) );

    return \@barcodes;
}

=head2 add_assay_failures

    my $barcodes2 = add_assay_failures( $failure_rate, $barcodes );

=cut

sub add_assay_failures {
    my ( $failure_rate, $barcodes ) =
      validate_pos( @_, $VAL_PROB, $VAL_BARCODES );
    return [
        map {
            [ map { rand() < $failure_rate ? 'X' : $_ } @$_ ]
        } @$barcodes
    ];
}

=head2 add_assay_failures_numeric

    my $numerics2 = add_assay_failures_numeric( $failure_rate, $numerics );

=cut

sub add_assay_failures_numeric {
    my ( $failure_rate, $numerics ) =
      validate_pos( @_, $VAL_PROB, $VAL_NUMERICS );
    return [
        map {
            [ map { rand() < $failure_rate ? 3 : $_ } @$_ ]
        } @$numerics
    ];
}

1;
