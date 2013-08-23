package COIL::Likelihood;

use strict;
use warnings;

use List::Util qw/ max sum /;
use List::MoreUtils qw/ pairwise /;

=head1 NAME

=cut

=head1 SUBROUTINES

=cut

=head2 posteriors

    my $posteriors = posteriors( \@likelihoods ); 

=cut

sub posteriors {
    my ($likelihoods) = @_;
    [ map { _posteriors($_) } @$likelihoods ];
}

sub _posteriors {
    my ($likelihoods) = @_;
    my $max = max @$likelihoods;
    my @unscaled = map { exp( $_ - $max ) } @$likelihoods;
    my $sum = sum @unscaled;
    return [ map { $_ / $sum } @unscaled ];
}

=head2 MAPs

    my $MAPs = MAPs( \@posteriors ); 

=cut

sub MAPs {
    return [
        map {
            my $d = $_;
            my $i = 0;
            $d->[$i] > $d->[$_] or $i = $_ for 1 .. $#$d;
            $i;
        } @{ $_[0] }
    ];
}

=head2 credible_intervals

    my $CIs = credible_intervals( \@posteriors, \@MAPs );
    my $CIs = credible_intervals( \@posteriors, \@MAPs, 0.95 );

=cut

sub credible_intervals {
    my ( $posteriors, $MAPs, $threshold ) = @_;
    $threshold ||= 0.95;

    return (
        [
            pairwise {
                my $lower = my $upper = $b;
                my $conf = $a->[$b];
                while ( $conf < $threshold ) {
                    if    ( $lower == 0 )    { $conf += $a->[ ++$upper ] }
                    elsif ( $upper == $#$a ) { $conf += $a->[ --$lower ] }
                    else {
                        if ( $a->[ $lower - 1 ] >= $a->[ $upper + 1 ] ) {
                            $conf += $a->[ --$lower ];
                        }
                        else { $conf += $a->[ ++$upper ] }
                    }
                }
                [ $lower, $upper, $conf ]
            }
            @$posteriors,
            @$MAPs
        ]
    );
}

1;
