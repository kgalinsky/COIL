package COIL::Tally::Pair;

use strict;
use warnings;

use Params::Validate;
use COIL::Validate ':val';

use COIL '_fh';

=head1 NAME

=head1 SYNOPSIS

=cut

=head1 CONSTRUCTORS

=cut

=head2 tally_numerics

    my $CTP = COIL::Tally::Pair->tally_numerics( \@numerics );

=cut

# A COIL::Tally::Pair object is a 2-dimensional array of counts
#
# $CTP->[$i][$j] = [ [ $NN_ij, $Nn_ij ], [ $nN_ij, $nn_ij ] ]

sub tally_numerics {
    my $class = shift;
    my ($numerics) = validate_pos( @_, $VAL_NUMERICS );

    my $self = bless [
        map {
            [ map { [ [ 0, 0 ], [ 0, 0 ] ] } @{ $numerics->[0] } ]
          } @{ $numerics->[0] }

    ], $class;

    foreach my $numeric (@$numerics) {
        next if ( grep { $_ == 2 } @$numeric );    # skip if barcode is poly

        # iterate through each SNP
        for ( my $i = 0 ; $i < @$numeric ; $i++ ) {
            my $ni = $numeric->[$i];
            next if ( $ni == 3 );                  # skip failed assays

            $self->[$i][$i][$ni][$ni]++;           # increment diagonal

            for ( my $j = $i + 1 ; $j < @$numeric ; $j++ ) {
                my $nj = $numeric->[$j];
                next if ( $nj == 3 );              # skip failed assays
                $self->[$i][$j][$ni][$nj]++;
            }
        }
    }

    # mirror
    for ( my $i = 1 ; $i < @$self ; $i++ ) {
        for ( my $j = 0 ; $j < $i ; $j++ ) {
            $self->[$i][$j] = $self->[$j][$i];
        }
    }

    return $self;
}

=head1 METHODS

=cut

=head2 Fisher

    my $p_values = $CTP->Fisher();

Perform Fisher's exact test on on each pair of SNPs. p-values are stored in a
2D symmetrical matrix with undef on the diagonal.

=cut

use Text::NSP::Measures::2D::Fisher::twotailed;

sub Fisher {
    my $self = shift;

    my @stats;
    for ( my $i = 0 ; $i < $#$self ; $i++ ) {
        for ( my $j = $i + 1 ; $j < @$self ; $j++ ) {
            my $n11 = $self->[$i][$j][0][0];
            my $n12 = $self->[$i][$j][0][1];
            my $n21 = $self->[$i][$j][1][0];
            my $n22 = $self->[$i][$j][1][1];

            my $n1p = $n11 + $n12;
            my $np1 = $n11 + $n21;
            my $npp = $n11 + $n12 + $n21 + $n22;

            $stats[$i][$j] = $stats[$j][$i] =
              Text::NSP::Measures::2D::Fisher::twotailed::calculateStatistic(
                n11 => $n11,
                n1p => $n1p,
                np1 => $np1,
                npp => $npp
              );
        }
    }

    return \@stats;
}

=head2 write

=cut

sub write {
    my $self = shift;
    my $fh = _fh( \@_, '>' );

    local $, = "\t";
    local $\ = "\n";
    for ( my $i = 0 ; $i < $#$self ; $i++ ) {
        for ( my $j = $i + 1 ; $j < @$self ; $j++ ) {
            my $t = $self->[$i][$j];
            print $i, $j,
                '['
              . $t->[0][0] . ':'
              . $t->[0][1] . '/'
              . $t->[1][0] . ':'
              . $t->[1][1] . ']';
        }
    }
}

1;
