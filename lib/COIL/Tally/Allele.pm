package COIL::Tally::Allele;

use strict;
use warnings;

use Carp;
use List::Util qw/ shuffle sum /;
use List::MoreUtils qw/ pairwise /;
use Math::Random qw/ random_beta /;

use Params::Validate;
use COIL::Validate qw/ :val :grab /;

=head1 NAME

	COIL::Tally::Allele - tally individual alleles and utilities

=head1 SYNOPSIS

Allele tallies are needed for converting barcodes to/from their numeric
representations as well as for computing likelihoods for each assay.

=cut

=head1 CONSTRUCTORS

=cut

=head2 new_from_barcodes

	my $CTA = COIL::Tally::Allele->new_from_barcodes($barcodes);

=cut

# A COIL::Tally::Allele object is an array of tuples:
#
# $CTA->[$i] = [ $A, $a, $N, $n ]

sub new_from_barcodes {
    my $class = shift;
    my ($barcodes) = validate_pos( @_, $VAL_BARCODES );

    my @seen;
    my @tallies;

    foreach my $barcode (@$barcodes) {
        my $seen_n;

        for ( my $i = 0 ; $i < @$barcode ; $i++ ) {
            my $snp = $barcode->[$i];
            if ( $snp =~ m/[ACGT]/ ) { $tallies[$i]{$snp} ||= 0 }
            elsif ( $snp eq 'N' ) { $seen_n = 1 }
        }

        unless ($seen_n) {
            for ( my $i = 0 ; $i < @$barcode ; $i++ ) {
                my $snp = $barcode->[$i];
                if ( $snp =~ m/[ACGT]/ ) { $tallies[$i]{$snp}++; }
            }
        }
    }

    return bless [ map { "${class}::Unit"->_new_from_hash($_) } @tallies ], $class;
}

=head2 new_rbeta

    my $CTA = COIL::Tally::Allele->new_rbeta( $n );
    my $CTA = COIL::Tally::Allele->new_rbeta( $n, $alpha, $beta );

Create a random tally with minor allele frequencies sampled from a Beta
distribution.

=cut

sub new_rbeta {
    my $class = shift;
    my ( $n, $a, $b ) = validate_pos(
        @_,
        { regex => qr/^[1-9]\d*$/ },
        (
            {
                default => 1,
                %$VAL_POS_REAL
            }
        ) x 2
    );

    return bless [
        map {
            my ( $A, $a ) = shuffle qw/ A C G T /;
            "${class}::Unit"->_new_from_p( $_, $A, $a );
        } random_beta( $n, $a, $b )
    ], $class;
}

=head2 new_rep

    my $CTA = COIL::Tally::Allele->new_rbeta( $n, $p );
    my $CTA = COIL::Tally::Allele->new_rbeta( $n, $p, $A, $a );

=cut

sub new_rep {
    my $class = shift;
    my ( $n, $p, $A, $a ) = validate_pos(
        @_, $VAL_POS_INT, $VAL_PROB,
        { %$VAL_SALLELE, default => 'A' },
        { %$VAL_SALLELE, default => 'T' },
    );

    return bless [ map { "${class}::Unit"->_new_from_p( $p, $A, $a ) } (0) x $n ],
      $class;
}

=head1 METHODS

=cut

=head2 barcodes2numerics

	my $numerics = $CTA->barcodes2numerics( \@barcodes );

=cut

sub barcodes2numerics {
    my $self = shift;
    my ($barcodes) = validate_pos( @_, $VAL_BARCODES );

    unless ( @$self == @{ $barcodes->[0] } ) {
        croak q{Tally length doesn't correspond to barcode length.};
    }

    my $struct = _barcode2numeric_struct( [ map { $_->[0] } @$self ] );
    return ( [ map { _barcode2numeric( $struct, $_ ) } @$barcodes ] );
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
    no warnings;
    [ pairwise { $a->[ ord($b) ] } @{ $_[0] }, @{ $_[1] } ];
}

=head2 numerics2barcodes

	my $barcodes = $CTA->numerics2barcodes( \@numerics );

=cut

sub numerics2barcodes {
    my $self = shift;
    my ($numerics) = validate_pos( @_, $VAL_NUMERICS );

    unless ( @$self == @{ $numerics->[0] } ) {
        croak q{Tally length doesn't correspond to barcode length.};
    }

    my $struct = _numeric2barcode_struct($self);
    return ( [ map { _numeric2barcode( $struct, $_ ) } @$numerics ] );
}

sub _numeric2barcode_struct {
    [ map { [ $_->[0], $_->[1], 'N', 'X' ] } @{ $_[0] } ];
}

sub _numeric2barcode {
    no warnings;
    [ pairwise { $a->[$b] } @{ $_[0] }, @{ $_[1] } ];
}

=head2 random_strain

    my $strain = $CTA->random_strain();

=cut

sub random_strain {
    [ map { $_->[ rand( $_->[2] + $_->[3] < $_->[3] ) ] } @{ $_[0] } ];
}

=head1 IO

=cut

=head2 read

    my $CTA = COIL::Tally::Allele->read( $fh );

=cut

sub read {
    my $class = shift;
    my $fh = grab_fh( @_, '<' );

    my $self = bless [], $class;

    local $/ = "\n";
    while ( local $_ = <$fh> ) {
        push @$self, [split];    # TODO validate tally lines
    }

    return $self;
}

=head2 write

    $CTA->write( $fh );

A direct output of the tally object.

    A   C   20  5
    G   T   34  6
    ...

=cut

sub write {
    my $self = shift;
    my $fh = grab_fh( \@_, '>' );

    local $, = "\t";
    local $\ = "\n";
    foreach my $allele (@$self) {
        print $fh @$allele;
    }
}

=head2 write_density

    $cta->write_density( $fh );

Write just the alleles and the density of the first allele. This will generally
be the major allele frequency.

    A   C   .8
    G   T   .85

=cut

sub write_density {
    my $self = shift;
    my $fh = grab_fh( \@_, '>' );

    local $, = "\t";
    local $\ = "\n";
    foreach my $allele (@$self) {
        print $fh @$allele[ 0, 1 ], $allele->[2] / sum( @$allele[ 2, 3 ] );
    }
}

package COIL::Tally::Allele::Unit;

use strict;
use warnings;

use Carp;

sub new {
    my $class = shift;
    bless [@_], $class;
}

sub _new_from_hash {
    my $class = shift;
    my ($tally) = @_;

    my @alleles = keys %$tally;

    if ( @alleles == 1 ) {
        carp qq{Monoallelic};
        $alleles[1] = '?';
    }

    elsif ( @alleles == 0 ) { croak qq{No alleles} }
    elsif ( @alleles > 2 ) {
        croak qq{Multiallelic not supported};
    }

    my ( $a, $b ) = @alleles;
    my ( $m, $n ) = @$tally{@alleles};
    $m ||= 0;
    $n ||= 0;

    my $self = $m > $n ? [ $a, $b, $m, $n ] : [ $b, $a, $n, $m ];
    bless $self, $class;
}

sub _new_from_p {
    my $class = shift;
    my ( $p, $A, $a ) = @_;

    my $N = 100 * $p;
    my $n = 100 - $N;

    my $self = $p >= .5 ? [ $A, $a, $N, $n ] : [ $a, $A, $n, $N ];
    bless $self, $class;
}

sub ref_allele  { $_[0][0] }
sub alt_allele  { $_[0][1] }
sub ref_count   { $_[0][2] }
sub alt_count   { $_[0][3] }
sub total_count { $_[0][2] + $_[0][3] }

sub p {
    ( $_[0]->ref_count + $_[1] ) / ( $_[0]->total_count + 2 * $_[1] );
}

1;
