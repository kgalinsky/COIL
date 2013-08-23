package COIL::Likelihood::Allele;

use strict;
use warnings;

use Params::Validate;

=head1 NAME

	COIL::Likelihood::Allele - allelic likelihoods

=head1 SYNOPSIS

Compute likelihoods for a barcode assuming independent alleles.

=cut

=head1 CONSTRUCTORS

=cut

=head2 tally2likelihood

	my $likelihood = COIL::Likelihood::Allele->tally2likelihood( $tally );
	my $likelihood = COIL::Likelihood::Allele->tally2likelihood(
	   $tally,
	   {
	       max_COI => 5,
	       padding => 0.5
	   }
	);

=cut

sub tally2likelihood {
    my $class = shift;
    my ( $tally, @p ) = validate_pos( @_, 1, { default => {} } );
    my %p = validate( @p, { max_COI => 5, padding => 0.5 } );

    my $max_COI = $p{max_COI};
    my $padding = $p{padding};

    my $self = bless [], $class;
    for ( my $i = 0 ; $i < $tally ; $i++ ) {
        my ( $p, $q ) = @{ $tally->[$i] }[ 2, 3 ];
        $p += $padding;
        $q += $padding;
        my $t = $p + $q;
        $p /= $t;
        $q /= $t;

        my $cum_log_p = my $log_p = log($p);
        my $cum_log_q = my $log_q = log($q);

        $self->[0][$i] = [ $log_p, $log_q, '-inf', 0 ];

        for ( my $c = 1 ; $c < $max_COI ; $c++ ) {
            $cum_log_p += $log_p;
            $cum_log_q += $log_q;
            my $log_n = log( 1 - exp($cum_log_p) - exp($cum_log_q) );
            $self->[$c][$i] = [ $cum_log_p, $cum_log_q, $log_n, 0 ];
        }
    }

    return ($self);
}

1;
