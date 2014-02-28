package COIL::Singular::Base;

use strict;
use warnings;

use Params::Validate;
use COIL::Validate ':val';

=head1 NAME

COIL::Barcode - interact with barcodes

=head1 SYNOPSIS

    my $barcodes = COIL::Barcodes->read( $filename_or_filehandle );
    my ( $barcodes, $lines ) = COIL::Barcodes->read( $filename_or_filehandle );
    my $numerics = $barcodes->to_numeric( \@major_alleles );
    foreach my $numeric ( @$numerics ) { print "$numeric\n" }

=head1 DESCRIPTION

This package contains 4 modules. Two of these are "singular" modules which
store one barcode (COIL::Barcode and COIL::Numeric) and two are "plural"
modules which store multiple barcodes (COIL::Barcodes and COIL::Numerics). The
plural modules are basically wrappers around the singular ones.

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

# Implementation in modules

=head2 read

    my $barcodes = COIL::Barcodes->read($filename);
    my $barcodes = COIL::Barcodes->read($filehandle);
    my ( $barcodes, $lines ) = COIL::Barcodes->read($filename);
    
    my $numerics = COIL::Numerics->read( $filename_or_filehandle );

Read barcodes from file. Creates an object of the plural class, which is a
container. If called in array context, it will also return the lines of the
file.

=cut

# Implemented in plural

=head1 METHODS

=cut

=head2 num_het

=head2 num_fail

=head2 eff_len

=head2 stats

=cut

sub eff_len {
    return scalar( @{ $_[0] } ) - $_[0]->num_fail;
}

sub stats {
    my $fail = $_[0]->num_fail;
    return [ $_[0]->num_het, $fail, scalar( @{ $_[0] } ) - $fail ]
}

=head2 to_string

    print $barcode->to_string();
    print $barcode;

    print $numeric->to_string();
    print $numeric;

Convert a barcode or numeric to a string.

=cut

# Generic implementation
use overload '""' => \&to_string;
sub to_string { join '', @{ $_[0] } }

=head2 add_assay_failures

    my $barcode2 = $barcode->add_assay_failures( $failure_rate );
    my $numeric2 = $numeric->add_assay_failures( $failure_rate );

    my $barcodes2 = $barcodes->add_assay_failures( $failure_rate );
    my $numerics2 = $numerics->add_assay_failures( $failure_rate );

Simulate assay failures.

=cut

# Generic implementation
sub add_assay_failures {
    shift->_add_assay_failures( validate_pos( @_, $VAL_PROB ) );
}

=head2 to_numeric / to_barcode

    my $numeric = $barcode->to_numeric( \@major alleles );
    my $barcode = $numeric->to_barcode( \@alleles );

    my $numerics = $barcodes->to_numeric( \@major_alleles );
    my $barcodes = $numerics->to_barcode( \@major_alleles );

=cut

# Implementation split between other base and derived classes

package COIL::Plural::Base;    # plural base class

use strict;
use warnings;

use Params::Validate;
use COIL::Validate qw/ :val :grab /;

# Generic constructors
sub new {
    my $class = shift;
    ( my $singular = $class ) =~ s/s$//;
    my ($self) =
      validate_pos( @_, { type => Params::Validate::ARRAYREF } );

    foreach my $el (@$self) {
        $el = $singular->new($el) unless ( ref($el) eq $singular );
    }

    bless $self, $class;
}

sub read {
    my $class = shift;
    ( my $singular = $class ) =~ s/s$//;
    my $lines = grab_lines( \@_ );
    my $self =
      $class->new( [ map { $singular->new_str( (split)[-1] ) } @$lines ] );

    return wantarray ? ( $self, $lines ) : $self;
}

sub write {
    my $self = shift;
    my $fh = grab_fh( \@_, '>' );

    local $\ = "\n";
    foreach my $el (@$self) { print $el }
}

# wrapper
sub add_assay_failures {
    my $self = shift;
    my ($failure_rate) = validate_pos( @_, $VAL_PROB );
    ref($self)
      ->new( [ map { $_->_add_assay_failures($failure_rate) } @$self ] );
}

package COIL::Barcode::Base;    # base class for barcode classes

# COIL::Barcode::Base and COIL::Numeric::Base were created to separate
# to_numeric/to_barcode. Thus, COIL::Barcode wouldn't have a defunct to_barcode
# function. These functions were also shared between the singular and plural
# versions of the classes, which is why a base class was created.

use strict;
use warnings;

use Params::Validate;
use COIL::Validate ':val';

sub to_numeric {
    my $self = shift;
    my ($major_alleles) = validate_pos( @_, $VAL_STRAIN );
    $self->_to_numeric( _to_numeric_struct($major_alleles) );
}

sub _to_numeric_struct {
    [
        map {
            my @a;
            @a[ ( 65, 67, 71, 84, 78, 88 ) ] = ( (1) x 4, 2, 3 );
            $a[ ord($_) ] = 0;
            \@a;
        } @{ $_[0] }
    ];
}

package COIL::Numeric::Base;    # base class for numerics

use strict;
use warnings;

use Params::Validate;
use COIL::Validate ':val';

#<<< define num_het and num_fail 
sub num_het  { return scalar grep { $_ == 2 } @{ $_[0] } }
sub num_fail { return scalar grep { $_ == 3 } @{ $_[0] } }
#>>>

sub to_barcode {
    my $self = shift;

    my ($alleles) = validate_pos( @_, { type => Params::Validate::ARRAYREF } );
    foreach my $pair (@$alleles) {
        validate_pos( @$pair, (@$VAL_SALLELE) x 2, (0) x 2 );
    }

    $self->_to_barcode( _to_barcode_struct($alleles) );
}

sub _to_barcode_struct {
    [ map { [ $_->[0], $_->[1], 'N', 'X' ] } @{ $_[0] } ];
}

package COIL::Barcode;    # barcode class

use strict;
use warnings;

our @ISA = qw/ COIL::Singular::Base COIL::Barcode::Base /;

use Carp;

use Params::Validate;
use COIL::Validate ':all';

use List::MoreUtils 'pairwise';

# CONSTRUCTORS
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

# METHODS

#<<< define num_het and num_fail 
sub num_het  { return scalar grep { $_ eq 'N' } @{ $_[0] } }
sub num_fail { return scalar grep { $_ eq 'X' } @{ $_[0] } }
#>>>

sub _add_assay_failures {
    my ( $self, $failure_rate ) = @_;
    ref($self)->new( [ map { rand() < $failure_rate ? 'X' : $_ } @$self ] );
}

sub _to_numeric {
    no warnings 'once';
    COIL::Numeric->new(
        [ pairwise { $b->[ ord($a) ] } @{ $_[0] }, @{ $_[1] } ] );
}

package COIL::Numeric;    # numeric class

use strict;
use warnings;

our @ISA = qw/ COIL::Singular::Base COIL::Numeric::Base /;

use Params::Validate;
use COIL::Validate ':val';

use List::MoreUtils 'pairwise';

# CONSTRUCTORS
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

# METHODS
sub _add_assay_failures {
    my ( $self, $failure_rate ) = @_;
    ref($self)->new( [ map { rand() < $failure_rate ? 3 : $_ } @$self ] );
}

sub _to_barcode {
    no warnings 'once';
    COIL::Barocde->new( [ pairwise { $b->[$a] } @{ $_[0] }, @{ $_[1] } ] );
}

package COIL::Barcodes;    # barcodes class

use strict;
use warnings;

our @ISA = qw/ COIL::Plural::Base COIL::Barcode::Base /;

sub _to_numeric {
    COIL::Numerics->new( [ map { $_->_to_numeric( $_[1] ) } @{ $_[0] } ] );
}

package COIL::Numerics;    # numerics class

use strict;
use warnings;

our @ISA = qw/ COIL::Plural::Base COIL::Numeric::Base /;

sub _to_barcode {
    COIL::Barcodes->new( [ map { $_->_to_barcode( $_[1] ) } @{ $_[0] } ] );
}

1;
