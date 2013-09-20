package COIL::Likelihood::Pair;

use strict;
use warnings;

use Params::Validate;
use COIL::Validate ':val';

=head1 NAME

COIL::Likelihood::Pair

=head1 SYNOPSIS

=cut

=head1 CONSTRUCTORS

=cut

=head2 new_from_tally

=cut

# A COIL::Likelihood::Pair object has the following structure:
#
# $CLA->[$c][$i][$j][$gi][$gj] = log P(G_i=gi & G_j=gj|C=c+1)

sub new_from_tally {
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
      "${class}::Level"->_new_from_Ps( [ map { $_->P($padding) } @$CTP ] );
    push @$self, $L;

    for ( my $i = 1 ; $i < $max_COI ; $i++ ) {
        $L = $L->_increment($L1);
        push @$self, $L;
    }

    return $self;
}

=head2 add_error

=cut

sub add_error {
    my $self = shift;
    my ($e) = validate_pos(
        @_,
        {
            default => 0.05,
            %$VAL_PROB
        }
    );

    # $E->[$i][$j] = P(G*=i|G=j)
    my $e1 = 1 - $e;
    my $e2 = $e / 2;
    my @E  = ( [ $e1, $e2, $e2 ], [ $e2, $e1, $e2 ], [ $e2, $e2, $e1 ] );

    return bless [ map { $_->_add_error( \@E ) } @$self ], ref($self);
}

=head2 numeric_likelihood

    my $likelihood  = $CLA->numeric_likelihood( $numeric );

=head2 numerics_likelihoods

    my $likelihoods = $CLA->numerics_likelihoods( \@numerics );

Compute log likelihood for a numeric or an array of numerics.

=cut

sub numeric_likelihood {
    shift->_numeric_likelihood( validate_pos( @_, $VAL_NUMERIC ) );
}

sub numerics_likelihoods {
    my $self = shift;
    my ($numerics) = validate_pos( @_, $VAL_NUMERICS );

    return ( [ map { $self->_numeric_likelihood($_) } @$numerics ] );
}

sub _numeric_likelihood {
    [ map { $_->_numeric_likelihood( $_[1] ) } @{ $_[0] } ];
}

package COIL::Likelihood::Pair::Level;

use strict;
use warnings;

use base 'COIL::Pair';

use List::MoreUtils 'pairwise';

sub _new_from_Ps {
    my $class = shift;
    my ($Ps) = @_;
    bless [ map { COIL::Likelihood::Pair::Unit->_new_from_P($_) } @$Ps ],
      $class;
}

sub _increment {
    my ( $self, $level0 ) = @_;
    no warnings 'once';
    bless [ pairwise { $a->_increment($b) } @$self, @$level0 ], ref($self);
}

sub _numeric_likelihood {
    my ( $self, $numeric ) = @_;

    my $l = 0;

    my $k = 0;
    for ( my $i = 1 ; $i < @$numeric ; $i++ ) {
        my $ni = $numeric->[$i];
        for ( my $j = 0 ; $j < $i ; $j++ ) {
            my $nj = $numeric->[$j];
            $l += $self->[ $k++ ][$ni][$nj];
        }
    }

    return $l / $#$numeric;
}

# Propogate error function
sub _add_error {
    my ( $self, $E ) = @_;
    bless [ map { $_->_add_error($E) } @$self ], ref($self);
}

package COIL::Likelihood::Pair::Unit;

use strict;
use warnings;

sub _new_from_P {
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
            $P->[$a1][$a2] = mylog( $P->[$a1][$a2] );
        }
    }

    # hets
    foreach my $a ( 0, 1, 3 ) {
        $P->[$a][2] = $P->[2][$a] = '-inf';
    }
    $P->[2][2] = '-inf';

    bless $P, $class;
}

sub _increment {
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
          mylog(
            exp( $unit->[$a][3] ) -
              exp( $unit->[$a][0] ) -
              exp( $unit->[$a][1] ) );
        $unit->[2][$a] =
          mylog(
            exp( $unit->[3][$a] ) -
              exp( $unit->[0][$a] ) -
              exp( $unit->[1][$a] ) );
    }

    return $unit;
}

sub _add_error {
    my ( $self, $E ) = @_;

    my $unit = bless [], ref($self);
    foreach my $i ( 0 .. 2 ) {
        foreach my $j ( 0 .. 2 ) {
            my $l = 0;

            foreach my $a ( 0 .. 2 ) {
                foreach my $b ( 0 .. 2 ) {
                    $l += exp( $self->[$a][$b] ) * $E->[$a][$i] * $E->[$b][$j];
                }
            }

            $unit->[$i][$j] = mylog($l);
        }

        $unit->[$i][3] =
          mylog(
            exp( $self->[0][3] ) * $E->[0][$i] +
              exp( $self->[1][3] ) * $E->[1][$i] +
              exp( $self->[2][3] ) * $E->[2][$i] );

        $unit->[3][$i] =
          mylog(
            exp( $self->[3][0] ) * $E->[0][$i] +
              exp( $self->[3][1] ) * $E->[1][$i] +
              exp( $self->[3][2] ) * $E->[2][$i] );

    }

    $unit->[3][3] = 0;

    return $unit;
}

sub mylog { $_[0] == 0 ? '-inf' : log $_[0] }

use overload '""' => \&to_string;

sub to_string {
    '[' . join(
        '/',
        map {
            join( ':', @$_ )
        } @{ $_[0] }
    ) . ']';
}

1;
