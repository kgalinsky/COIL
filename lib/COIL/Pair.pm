package COIL::Pair;

use strict;
use warnings;

=head1 NAME

    COIL::Pair - base class for objects that deal with pairs of alleles

=head1 SYNOPSIS

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
