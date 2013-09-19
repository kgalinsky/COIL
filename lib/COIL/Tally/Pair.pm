package COIL::Tally::Pair;

use strict;
use warnings;

use base 'COIL::Pair';

use Params::Validate;
use COIL::Validate ':val';

=head1 NAME

=head1 SYNOPSIS

=cut

=head1 CONSTRUCTORS

=cut

=head2 tally_numerics

    my $CTP = COIL::Tally::Pair->tally_numerics( \@numerics );

=cut

# $CTP->[$k] = [ [ $NN_ij, $Nn_ij ], [ $nN_ij, $nn_ij ] ]

sub tally_numerics {
    my $class = shift;
    my ($numerics) = validate_pos( @_, $VAL_NUMERICS );

    my $n = $#{ $numerics->[0] };
    my $self =
      bless [ map { "${class}::Unit"->new } (0) x ( $n * ( $n + 1 ) / 2 ) ],
      $class;

    foreach my $numeric (@$numerics) {
        next if ( grep { $_ == 2 } @$numeric );    # skip if barcode is poly

        # iterate through each SNP
        for ( my $i = 1 ; $i < @$numeric ; $i++ ) {
            my $ni = $numeric->[$i];
            next if ( $ni == 3 );                  # skip failed assays

            my $offset = $i * ( $i - 1 ) / 2;

            for ( my $j = 0 ; $j < $i ; $j++ ) {
                my $nj = $numeric->[$j];
                next if ( $nj == 3 );              # skip failed assays
                $self->[ $offset + $j ][$ni][$nj]++;
            }
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
    return bless [
        map {
            my $n11 = $_->[0][0];
            my $n12 = $_->[0][1];
            my $n21 = $_->[1][0];
            my $n22 = $_->[1][1];

            my $n1p = $n11 + $n12;
            my $np1 = $n11 + $n21;
            my $npp = $n11 + $n12 + $n21 + $n22;

            Text::NSP::Measures::2D::Fisher::twotailed::calculateStatistic(
                n11 => $n11,
                n1p => $n1p,
                np1 => $np1,
                npp => $npp
            );
        } @$self
    ], 'COIL::Pair';
}

package COIL::Tally::Pair::Unit;

use strict;
use warnings;

sub new {
    my $class = shift;
    bless [ [ 0, 0 ], [ 0, 0 ] ], $class;
}

use overload '""' => \&to_string;
sub to_string { "[$_[0][0][0]:$_[0][0][1]/$_[0][1][0]:$_[0][1][1]]" }

1;
