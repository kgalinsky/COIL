package COIL::Likelihood::Pair;

use strict;
use warnings;

use base 'COIL::Pair';

use Params::Validate;
use COIL::Validate ':val';

=head1 NAME

COIL::Likelihood::Pair

=head1 SYNOPSIS

=cut

=head1 CONSTRUCTORS

=cut

=head2 tally2likelihood

=cut

# A COIL::Likelihood::Pair object has the following structure:
#
# $CLA->[$c][$i][$j][$gi][$gj] = log P(G_i=gi & G_j=gj|C=c+1)

sub tally2likelihood {
    my $class = shift;
    my ( $CTP, @p ) = validate_pos( @_, 1, { default => {} } );
    my %p = validate(
        @p,
        {
            max_COI => { default => 5,   %$VAL_POS_INT },
            padding => { default => 0.5, %$VAL_NON_NEG_REAL }
        }
    );

    my $max_COI = $p{max_COI};
    my $padding = $p{padding};

    my $self = bless [], $class;
    my $L1 = my $L =
      "${class}::Level"->new_Ps( [ map { $_->density($padding) } @$CTP ] );
      push @$self, $L;

    for (my $i = 1; $i < $max_COI; $i++) {
        $L = $L->increment($L1);
        push @$self, $L;
    }
    
    return $self;
}

package COIL::Likelihood::Pair::Level;

use strict;
use warnings;

use List::MoreUtils 'pairwise';

sub new_Ps {
    my $class = shift;
    my ($Ps) = @_;
    bless [ map { COIL::Likelihood::Pair::Unit->new_P($_) } @$Ps ], $class;
}

sub increment {
    my ( $self, $level0 ) = @_;
    no warnings 'once';
    bless [ pairwise { $a->increment($b) } @$self, @$level0 ], ref($self);
}

package COIL::Likelihood::Pair::Unit;

use strict;
use warnings;

sub new_P {
    my $class = shift;
    my ($P) = @_;

    # marginals
    foreach my $a ( 0, 1 ) {
        $P->[$a][3] = $P->[$a][0] + $P->[$a][1];
        $P->[3][$a] = $P->[0][$a] + $P->[1][$a];
    }
    $P->[3][3] = 1;

    # log
    foreach my $a1 ( 0, 1, 3 ) {
        foreach my $a2 ( 0, 1, 3 ) {
            $P->[$a1][$a2] = log $P->[$a1][$a2];
        }
    }

    # hets
    foreach my $a ( 0, 1, 3 ) {
        $P->[$a][2] = $P->[2][$a] = '-inf';
    }

    bless $P, $class;
}

sub increment {
    my ( $self, $unit0 ) = @_;
    my $unit = bless [], ref($self);

    # P(Gi=0,Gj=0|C=c) = P(Gi=0,Gj=0|C=1)^c
    foreach my $a1 ( 0, 1, 3 ) {
        foreach my $a2 ( 0, 1, 3 ) {
            $unit->[$a1][$a2] =
              $self->[$a1][$a2] + $unit0->[$a1][$a2];
        }
    }

    # P(Gi=0,Gj=2|C=c) = P(Gi=0,Gj=3) - P(Gi=0,Gj=0) - P(Gi=0,Gj=1)
    foreach my $a ( 0, 1, 3, 2 ) {
        $unit->[$a][2] =
          log(
            exp( $unit->[$a][3] ) -
              exp( $unit->[$a][0] ) -
              exp( $unit->[$a][1] ) );
        $unit->[2][$a] =
          log(
            exp( $unit->[3][$a] ) -
              exp( $unit->[0][$a] ) -
              exp( $unit->[1][$a] ) );
    }

    return $unit;
}

1;
