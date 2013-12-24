package COIL::Tally::Pair;

use strict;
use warnings;

use parent 'COIL::Pair';

use Params::Validate;
use COIL::Validate ':val';

=head1 NAME

COIL::Tally::Pair - tally pairs of alleles in barcodes

=head1 SYNOPSIS

    my $tally    = COIL::Tally::Pair->new_from_numerics( \@numerics );
    my $p_values = $tally->fisher();

=cut

=head1 CONSTRUCTORS

=cut

=head2 new_from_numerics

    my $tally = COIL::Tally::Pair->new_from_numerics( \@numerics );

=cut

sub new_from_numerics {
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

=head2 fisher

    my $p_values = $tally->fisher();

Perform Fisher's exact test on on each pair of SNPs. They are returned in the
same order as this object as well as what is specified in COIL::Pair.

=cut

sub fisher {
    return bless [ map { $_->fisher() } @{ $_[0] } ], 'COIL::Pair';
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

sub refref_count { $_[0][0][0] }
sub refalt_count { $_[0][0][1] }
sub altref_count { $_[0][1][0] }
sub altalt_count { $_[0][1][1] }
sub total_count  { $_[0][0][0] + $_[0][0][1] + $_[0][1][0] + $_[0][1][1] }

sub P {
    my $total = $_[0]->total_count + 4 * $_[1];

    return [
        [ ( $_[0][0][0] + $_[1] ) / $total, ( $_[0][0][1] + $_[1] ) / $total ],
        [ ( $_[0][1][0] + $_[1] ) / $total, ( $_[0][1][1] + $_[1] ) / $total ]
    ];
}

use Text::NSP::Measures::2D::Fisher::twotailed;

sub fisher {
    my $n11 = $_[0]->[0][0];
    my $n12 = $_[0]->[0][1];
    my $n21 = $_[0]->[1][0];
    my $n22 = $_[0]->[1][1];

    my $n1p = $n11 + $n12;
    my $np1 = $n11 + $n21;
    my $npp = $n11 + $n12 + $n21 + $n22;

    Text::NSP::Measures::2D::Fisher::twotailed::calculateStatistic(
        n11 => $n11,
        n1p => $n1p,
        np1 => $np1,
        npp => $npp
    );
}

1;
