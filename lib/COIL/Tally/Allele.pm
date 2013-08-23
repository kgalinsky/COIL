package COIL::Tally::Allele;

use strict;
use warnings;

use Carp;
use List::MoreUtils 'pairwise';

use Params::Validate;
use COIL::Validate ':val';

=head1 NAME

	COIL::Tally::Allele - tally individual alleles and utilities

=head1 SYNOPSIS

Allele tallies are needed for converting barcodes to/from their numeric
representations as well as for computing likelihoods for each assay.

=cut

=head1 CONSTRUCTORS

=cut

=head2 tally_barcodes

	my $tally = COIL::Tally::Allele->tally_barcodes($barcodes);

=cut

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
        croak qq{Site "$i" not biallelic} unless ( @alleles == 2 );

        my ( $a, $b ) = @alleles;
        my ( $m, $n ) = @{ $tallies[$i] }{@alleles};
        $m ||= 0;
        $n ||= 0;

        push @$self, $m > $n ? [ $a, $b, $m, $n ] : [ $b, $a, $n, $m ];
    }

    return ($self);
}

=head1 METHODS

=cut

=head2 barcodes2numerics

	my $numerics = $tally->barcodes2numerics( \@barcodes );

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

1;
