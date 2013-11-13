package COIL::Barcode;

use strict;
use warnings;

use Carp;

use Params::Validate;
use COIL::Validate ':all';

=head1 CONSTRUCTORS

=cut

=head2 new

    my $barcode = COIL::Barcode->new( \@barcode );

=head2 new_str

    my $barcode = COIL::Barcode->new_str( $barcode_string );

=cut

sub new {
    my $class = shift;
    my ($barcode) = validate_pos(@_, $VAL_BARCODE);
    bless $barcode, $class;
}

sub new_str {
    my $class = shift;
    my ($barcode_str)= validate_pos(@_, $VAL_BARCODE_STR);
    bless [ split m//, $barcode_str], $class;
}

=head1 FUNCTIONS

=head2 read_barcodes

    my $barcodes = read_barcodes($filename);
    my $barcodes = read_barcodes($filehandle);

Read barcodes from file.

=cut

sub read_barcodes {
    my $lines = grab_lines( \@_ );
    my @barcodes = 
      map  { [ split m// ] }    # split into array of chars
      map  { (split)[-1] }      # grab last column
      @$lines;

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

package COIL::Barcodes;

use strict;
use warnings;

package COIL::Numeric;



1;
