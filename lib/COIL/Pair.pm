package COIL::Pair;

use strict;
use warnings;

use COIL '_fh';

=head1 NAME

    COIL::Pair - base class for objects that deal with pairs of alleles

=head1 SYNOPSIS

    $pair->write();

=head1 DESCRIPTION

It is sometimes useful to deal with pairs of alleles in COIL. There were two
options for handling this. One is to make a 2D array. The other possibility is
a 1D array with a special indexing. In order for map to work properly, the
latter is preferable. It has the added benefit of using less memory and being
a "simpler" data structure.

A COIL::Pair object is an array of counts stores the bottom-left off-diagonal
elements of a 2D matrix. The numbering is:

         0   1   2   3
      +---+---+---+---+---+
    0 | X |   |   |   |...|
      +---+---+---+---+---+
    1 | 0 | X |   +   +...|
      +---+---+---+---+---+
M = 2 | 1 | 2 | X |   |...|
      +---+---+---+---+---+
    3 | 3 | 4 | 5 | X |...|
      +---+---+---+---+---+
      |...|...|...|...|...|
      +---+---+---+---+---+

To map M[i][j] = P[k]:

    k = i(i-1)/2+j

=cut

=head2 write

    $pair->write();
    $pair->write( $filename );
    $pair->write( $filehandle );
    COIL::Pair::write( $pair_array );

=cut

sub write {
    my $self = shift;
    my $fh = _fh( \@_, '>' );

    local $, = "\t";
    local $\ = "\n";

    my ( $i, $j ) = ( 1, 0 );
    foreach my $t (@$self) {
        print $i, $j, $t;

        $j++;
        if ( $i == $j ) { $i++; $j = 0; }
    }
}

1;
