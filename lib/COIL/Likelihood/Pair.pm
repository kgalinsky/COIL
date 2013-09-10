package COIL::Likelihood::Pair;

use strict;
use warnings;

use Params::Validate;

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
            max_COI => { default => 5 },
            padding => { default => 0.5 }
        }
    );

    my $max_COI = $p{max_COI};
    my $padding = $p{padding};

    my $self = bless [], $class;
    for ( my $i = 1 ; $i < @$CTP ; $i++ ) {
        for ( my $j = 0 ; $j < $i ; $j++ ) {
            my @P;
            my $t = 0;

            # Initialize counts + padding and total
            foreach my $a1 ( 0, 1 ) {
                foreach my $a2 ( 0, 1 ) {
                    $t += $P[$a1][$a2] = $CTP->[$i][$j][$a1][$a2] + $padding;
                }
            }

            # Convert to probabilities
            foreach my $a1 ( 0, 1 ) {
                foreach my $a2 ( 0, 1 ) {
                    $P[$a1][$a2] /= $t;
                }
            }

            # Compute marginals
            foreach my $a ( 0, 1 ) {
                $P[$a][3] = $P[$a][0] + $P[$a][1];
                $P[3][$a] = $P[0][$a] + $P[1][$a];
            }
            $P[3][3] = 1;

            # Take the log
            foreach my $a1 ( 0, 1, 3 ) {
                foreach my $a2 ( 0, 1, 3 ) {
                    $P[$a1][$a2] = log $P[$a1][$a2];
                }
            }

            # Set the poly to -inf
            foreach my $a ( 0, 1, 3 ) {
                $P[$a][2] = $P[2][$a] = '-inf';
            }

            # Store
            $self->[0][$i][$j] = \@P;

            # Copy relevant cells to Q
            my @Q;
            foreach my $a1 ( 0, 1, 3 ) {
                foreach my $a2 ( 0, 1, 3 ) {
                    $Q[$a1][$a2] = $P[$a1][$a2];
                }
            }

            for ( my $c = 1 ; $c < $max_COI ; $c++ ) {

                # P(Gi=0,Gj=0|C=c) = P(Gi=0,Gj=0|C=1)^c
                foreach my $a1 ( 0, 1, 3 ) {
                    foreach my $a2 ( 0, 1, 3 ) {
                        $Q[$a1][$a2] += $P[$a1][$a2];
                    }
                }

                # P(Gi=0,Gj=2|C=c) = P(Gi=0,Gj=3) - P(Gi=0,Gj=0) - P(Gi=0,Gj=1)
                foreach my $a ( 0, 1, 3, 2 ) {
                    $Q[$a][2] =
                      log(
                        exp( $Q[$a][3] ) -
                          exp( $Q[$a][0] ) -
                          exp( $Q[$a][1] ) );
                    $Q[2][$a] =
                      log(
                        exp( $Q[3][$a] ) -
                          exp( $Q[0][$a] ) -
                          exp( $Q[1][$a] ) );
                }

                # Copy Q to R and store
                my @R;
                foreach my $a1 ( 0 .. 3 ) {
                    foreach my $a2 ( 0 .. 3 ) {
                        $R[$a1][$a2] = $Q[$a1][$a2];
                    }
                }
                $self->[$c][$i][$j] = \@R;
            }
        }
    }

    return $self;
}

1;
