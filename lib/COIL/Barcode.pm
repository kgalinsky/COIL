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
    my $numeric = COIL::Numeric->new( \@numeric );

Create a barcode or numeric from an arrayref of alleles or numeric alleles.

=head2 new_str

    my $barcode = COIL::Barcode->new_str( $barcode_string );
    my $numeric = COIL::Numeric->new_str( $numeric_string );

Create a barcode or numeric from a barcode or numeric string.

=cut

sub new {
    my $class = shift;
    my ($barcode) = validate_pos( @_, $VAL_BARCODE );
    bless $barcode, $class;
}

sub new_str {
    my $class = shift;
    my ($barcode_str) = validate_pos( @_, $VAL_BARCODE_STR );
    bless [ split m//, $barcode_str ], $class;
}

=head2 read

    my $barcodes = COIL::Barcode->read($filename);
    my $barcodes = COIL::Barcode->read($filehandle);
    
    my $numerics = COIL::Numeric->read( $filename_or_filehandle );

Read barcodes from file. Creates an object of the plural class, which is a
container.

=cut

sub read {
    my $class  = shift;
    my $plural = $class . 's';
    my $lines  = grab_lines( \@_ );
    $plural->new( [ map { $class->new_str( (split)[-1] ) } @$lines ] );
}

=head2 add_assay_failures

    my $barcode2 = $barcode->add_assay_failures( $failure_rate );
    my $numeric2 = $numeric->add_assay_failures( $failure_rate );

    my $barcodes2 = $barcodes->add_assay_failures( $failure_rate );
    my $numerics2 = $numerics->add_assay_failures( $failure_rate );

Simulate assay failures.

=cut

sub add_assay_failures {
    shift->_add_assay_failures( validate_pos( @_, $VAL_PROB ) );
}

sub _add_assay_failures {
    my ( $self, $failure_rate ) = @_;
    ref($self)->new( [ map { rand() < $failure_rate ? 'X' : $_ } @$self ] );
}

package COIL::Numeric;

use strict;
use warnings;

our @ISA = 'COIL::Barcode';

use Params::Validate;
use COIL::Validate ':val';

sub new {
    my $class = shift;
    my ($barcode) = validate_pos( @_, $VAL_NUMERIC );
    bless $barcode, $class;
}

sub new_str {
    my $class = shift;
    my ($barcode_str) = validate_pos( @_, $VAL_NUMERIC_STR );
    bless [ split m//, $barcode_str ], $class;
}

sub _add_assay_failures {
    my ( $self, $failure_rate ) = @_;
    ref($self)->new( [ map { rand() < $failure_rate ? 3 : $_ } @$self ] );
}

package COIL::Barcodes;

use strict;
use warnings;

use Params::Validate;
use COIL::Validate ':val';

sub new {
    my $class = shift;
    ( my $singular = $class ) =~ s/s$//;
    my ($self) = validate_pos( @_, { type => Params::Validate::ARRAYREF } );
    validate_pos( @$self, ( { isa => $singular } ) x @$self );
    bless $self, $class;
}

sub add_assay_failures {
    my $self = shift;
    my ($failure_rate) = validate_pos( @_, $VAL_PROB );
    ref($self)
      ->new( [ map { $_->_add_assay_failures($failure_rate) } @$self ] );
}

package COIL::Numerics;

use strict;
use warnings;

our @ISA = 'COIL::Barcodes';

1;
