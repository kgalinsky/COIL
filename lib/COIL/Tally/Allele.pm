package COIL::Tally::Allele;

use strict;
use warnings;

use Carp;
use Scalar::Util qw/ looks_like_number /;
use List::Util qw/ shuffle sum /;
use List::MoreUtils qw/ pairwise /;
use Math::Random qw/ random_beta /;

use Params::Validate;
use COIL::Validate ':val';

use COIL '_fh';

=head1 NAME

	COIL::Tally::Allele - tally individual alleles and utilities

=head1 SYNOPSIS

Allele tallies are needed for converting barcodes to/from their numeric
representations as well as for computing likelihoods for each assay.

=cut

=head1 CONSTRUCTORS

=cut

=head2 tally_barcodes

	my $CTA = COIL::Tally::Allele->tally_barcodes($barcodes);

=cut

# A COIL::Tally::Allele object is an array of tuples:
#
# $CTA->[$i] = [ $A, $a, $N, $n ]

sub tally_barcodes {
    my $class = shift;
    my ($barcodes) = validate_pos( @_, $VAL_BARCODES );

    my @seen;
    my @tallies;

    foreach my $barcode (@$barcodes) {
        my $seen_n;

        for ( my $i = 0 ; $i < @$barcode ; $i++ ) {
            my $snp = $barcode->[$i];
            if ( $snp =~ m/[ACGT]/ ) { $seen[$i]{$snp} = 1 }
            elsif ( $snp eq 'N' ) { $seen_n = 1 }
        }

        unless ($seen_n) {
            for ( my $i = 0 ; $i < @$barcode ; $i++ ) {
                my $snp = $barcode->[$i];
                if ( $snp =~ m/[ACGT]/ ) { $tallies[$i]{$snp}++; }
            }
        }
    }

    my $self = bless [], $class;

    for ( my $i = 0 ; $i < @seen ; $i++ ) {
        my @alleles = keys %{ $seen[$i] };
        if ( @alleles == 1 ) {
            carp qq{Site $i not biallelic};
            $alleles[1] = '?';
        }
        elsif ( @alleles == 0 ) { croak qq{No alleles at site $i} }
        elsif ( @alleles > 2 ) {
            croak qq{Multiallelic not supported at site $i};
        }

        my ( $a, $b ) = @alleles;
        my ( $m, $n ) = @{ $tallies[$i] }{@alleles};
        $m ||= 0;
        $n ||= 0;

        push @$self, $m > $n ? [ $a, $b, $m, $n ] : [ $b, $a, $n, $m ];
    }

    return ($self);
}

=head2 random_tally

    my $CTA = COIL::Tally::Allele->random_tally( $n );
    my $CTA = COIL::Tally::Allele->random_tally( $n, $alpha, $beta );

=cut

sub random_tally {
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
            $_ *= 100;
            my ( $N, $n ) = $_ > 50 ? ( $_, 100 - $_ ) : ( 100 - $_, $_ );
            my ( $A, $a ) = shuffle qw/ A C G T /;
            [ $A, $a, $N, $n ];
        } random_beta( $n, $a, $b )
    ], $class;
}

=head2 uniform_tally

    my $CTA = COIL::Tally::Allele->random_tally( $n, $p );
    my $CTA = COIL::Tally::Allele->random_tally( $n, $p, $A, $a );

=cut

sub uniform_tally {
    my $class = shift;
    my ( $n, $p, $A, $a ) = validate_pos(
        @_, $VAL_POS_INT, $VAL_PROB,
        { %$VAL_NUC, default => 'A' },
        { %$VAL_NUC, default => 'T' },
    );
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
    my $fh = _fh( @_, '<' );

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
    my $fh = _fh( \@_, '>' );

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
    my $fh = _fh( \@_, '>' );

    local $, = "\t";
    local $\ = "\n";
    foreach my $allele (@$self) {
        print $fh @$allele[ 0, 1 ], $allele->[2] / sum( @$allele[ 2, 3 ] );
    }
}

1;
