package COIL::Likelihood::Allele;

use strict;
use warnings;

use List::Util qw/ sum /;
use List::MoreUtils qw/ pairwise /;

use Params::Validate;
use COIL::Validate ':val';

use COIL '_fh';

=head1 NAME

	COIL::Likelihood::Allele - allelic likelihoods

=head1 SYNOPSIS

Compute likelihoods for a barcode assuming independent alleles.

=cut

=head1 CONSTRUCTORS

=cut

=head2 tally2likelihood

	my $CLA = COIL::Likelihood::Allele->tally2likelihood( $tally );
	my $CLA = COIL::Likelihood::Allele->tally2likelihood(
	   $tally,
	   {
	       max_COI => 5,
	       padding => 0.5
	   }
	);

=cut

# A COIL::Likelihood::Allele object has the following structure:
#
# $CLA->[$c][$i][$g] = log P(G_i=g|C=c+1)
#
# This was chosen over $CLA->[$i][$c][$g]. Grouping things by SNP first has
# benefits, but the key advantage of clustering by COI is that computing the
# genotype likelihood (log P(G|C=c)) can be done by zipping with $CLA->[$c].

sub tally2likelihood {
    my $class = shift;
    my ( $tally, @p ) = validate_pos( @_, 1, { default => {} } );
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
    for ( my $i = 0 ; $i < @$tally ; $i++ ) {
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

=head1 METHODS

=cut

=head2 add_error

    my $CLA_E = $CLA->add_error();
    my $CLA_E = $CLA->add_error( 0.05 );

=cut

sub add_error {
    my $self = shift;
    my ($e) = @_;

    $e ||= .05;

    # $E->[$i][$j] = P(G*=i|G=j)
    my $E = [
        [ 1 - 2 * $e, $e,         $e,         0 ],
        [ $e,         1 - 2 * $e, $e,         0 ],
        [ $e,         $e,         1 - 2 * $e, 0 ],
        [ 0,          0,          0,          1 ]
    ];

    my $CLA_E = bless [], ref($self);
    for ( my $c = 0 ; $c < @$self ; $c++ ) {
        for ( my $i = 0 ; $i < @{ $self->[$c] } ; $i++ ) {
            $CLA_E->[$c][$i] = [
                map {
                    log( sum pairwise { $a * exp($b) } @$_,
                        @{ $self->[$c][$i] } )
                } @$E
            ];
        }
    }

    return ($CLA_E);
}

=head2 numerics_likelihoods

    my $likelihoods = $CLA->numerics_likelihoods( $numerics );

=cut

sub numerics_likelihoods {
    my $self = shift;
    my ($numerics) = validate_pos( @_, $VAL_NUMERICS );

    return ( [ map { $self->_numeric_likelihoods($_) } @$numerics ] );
}

sub _numeric_likelihoods {
    my ( $self, $numeric ) = @_;

    no warnings;
    return (
        [
            map {
                sum pairwise { $a->[$b] } @$_, @$numeric
            } @$self
        ]
    );
}

=head2 write

    $CLA->write( $fh, $digits);

=cut

sub write {
    my $self     = shift;
    my $fh       = _fh( \@_, '>' );
    my ($digits) = @_;
    $digits ||= 2;

    for ( my $i = 0 ; $i < @{ $self->[0] } ; $i++ ) {
        for ( my $c = 0 ; $c < @$self ; $c++ ) {
            if ( $c != 0 ) { print $fh "\t" }
            print $fh join( '|',
                map { sprintf "\%.${digits}f", exp($_) }
                  @{ $self->[$c][$i] }[ 0 .. 2 ] );
        }
        print $fh "\n";
    }
}

1;
