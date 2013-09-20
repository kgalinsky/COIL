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

Compute log likelihoods for a barcode assuming independent alleles.

    my $CLA = COIL::Likelihood::Allele->tally2likelihood( $CTA );
    my $CLA_E = $CLA->add_error( $e );
    my $likelihoods = $CLA_E->numerics_likelihoods( $numerics );

=head1 DESCRIPTION

Given a tally object, create a likelihood object that can be used to compute
a likelihood at a range of COIs for a barcode. The likelihoods that are stored
and computed are actually log likelihoods because they are not as succeptible
to rounding errors.

The following equation is used for calculating the likelihood:

    L(C|G) = P(G|C) = prod P(G_i|C)
    l(C|G) = log(L(C|G)) = sum log(P(G_i|C))

G_i is a random variable that has the following values:

    0: Homozygous reference (let A be when G_i=0)
    1: Homozygous alternate (let a be when G_i=1)
    2: Heterozygous         (let N be when G_i=2)
    3: Failed assay         (let X be when G_i=3)

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

Create the likelihood calculation object. This is a log of the following
vector:

        / P(A|C) \
    L = | P(a|C) |
        \ P(N|C) /

There is a fourth element, log(P(X|C))=0 since failed assays do not tell us any
information about COI. To compute these elements, let B be a count of ref calls
in the different strains.

    B ~ Binomial(C, p), p = allele frequency of ref

    A = I(B=C)
    a = I(B=0)
    N = I(B>0, B<C)

    P(A) = P(B=C) = p^C
    P(a) = P(B=0) = (1-p)^C
    P(N) = 1 - P(A) - P(a)

Since we are doing everything on a log scale, the first two equations become:

    log P(A) = C * log(p)
    log P(a) = C * log(1-p)

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
            max_COI => { default => 5,   %$VAL_POS_INT },
            padding => { default => 0.5, %$VAL_NON_NEG_REAL }
        }
    );

    my $max_COI = $p{max_COI};
    my $padding = $p{padding};

    my @ladders =
      map { "${class}::Unit"->ladder( $_->p($padding), $max_COI ) } @$tally;
    return bless [
        map {
            my $c = $_;
            "${class}::Level"->new( [ map { $_->[$c] } @ladders ] )
        } ( 0 .. $max_COI - 1 )
      ],
      $class;
}

=head1 METHODS

=cut

=head2 add_error

    my $CLA_E = $CLA->add_error();
    my $CLA_E = $CLA->add_error( 0.05 );
    my $CLA_E = $CLA->add_error( $E );

This is computed by doing a matrix-vector multiplication. The matrix may be
specified directly (E) or using just a simple error rate (e). Let G_i* be the
observed allele and A*/a*/N* be the associated events from before. The matrix
is just a stochastic matrix.

    E[i][j] = P(G_i*=j|G_i=i)

        / P(A*|A), P(a*|A), P(N*|A) \
    E = | P(A*|a), P(a*|a), P(N*|a) |
        \ P(A*|N), P(a*|N), P(N*|N) /

        / 1-e, e/2, e/2 \
      = | e/2, 1-e, e/2 |, simple error rate case
        \ e/2, e/2, 1-e /

To compute the new likelihood (L*) given the old one (L), simply multiply:

    L* = LE

Note that in this representation, the likelihood vectors are 1x3 rather than
3x1. To use 3x1 vectors:

    L* = E'L

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
    my @ET = ( [ $e1, $e2, $e2 ], [ $e2, $e1, $e2 ], [ $e2, $e2, $e1 ] );

    return bless [ map { $_->add_error( \@ET ) } @$self ], ref($self);
}

=head2 numerics_likelihoods

    my $likelihoods = $CLA->numerics_likelihoods( $numerics );

Compute log likelihoods at different COIs for each numeric.

=cut

sub numerics_likelihoods {
    my $self = shift;
    my ($numerics) = validate_pos( @_, $VAL_NUMERICS );

    return ( [ map { $self->_numeric_likelihood($_) } @$numerics ] );
}

sub _numeric_likelihood {
    [ map { $_->numeric_likelihood( $_[1] ) } @{ $_[0] } ];
}

=head2 random_numeric

    my $numeric = $CLA->random_numeric( $COI );

Generate a random numeric at a particular COI.

=cut

sub random_numeric {
    my $self = shift;
    my ($COI) = validate_pos( @_, $VAL_POS_INT );
    return [ map { _random_numeric($_) } @{ $self->[ $COI - 1 ] } ];
}

sub _random_numeric {
    my ($l) = @_;
    my $r = rand(1);

    for ( my $n = 0 ; $n < 3 ; $n++ ) {
        my $p = exp( $l->[$n] );
        return $n if ( $r < $p );
        $r -= $p;
    }

    croak(
        sprintf(
            'Random value greater than probabilities: %g %g %g',
            $l->[0], $l->[1], $l->[2]
        )
    );
}

=head2 write

    $CLA->write( $fh, $digits);

=cut

sub write {
    my $self = shift;
    my $fh = _fh( \@_, '>' );
    my ($digits) = validate_pos( @_, { default => 2, %$VAL_POS_INT } );

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

package COIL::Likelihood::Allele::Level;

use strict;
use warnings;

use List::Util 'sum';
use List::MoreUtils 'pairwise';

sub new { bless $_[1], $_[0] }

sub numeric_likelihood {
    my ( $self, $numeric ) = @_;
    no warnings 'once';
    sum pairwise { $a->[$b] } @$self, @$numeric;
}

sub add_error {
    my ( $self, $ET ) = @_;
    bless [ map { $_->add_error($ET) } @$self ], ref($self);
}

package COIL::Likelihood::Allele::Unit;

use strict;
use warnings;

use List::Util 'sum';
use List::MoreUtils 'pairwise';

sub ladder {
    my $class = shift;
    my ( $p, $max_COI ) = @_;
    my $q = 1 - $p;

    my $cum_log_p = my $log_p = log($p);
    my $cum_log_q = my $log_q = log($q);

    my @selves;
    push @selves, bless [ $log_p, $log_q, '-inf', 0 ], $class;

    for ( my $c = 1 ; $c < $max_COI ; $c++ ) {
        $cum_log_p += $log_p;
        $cum_log_q += $log_q;
        my $log_n = log( 1 - exp($cum_log_p) - exp($cum_log_q) );
        push @selves, bless [ $cum_log_p, $cum_log_q, $log_n, 0 ], $class;
    }

    return \@selves;
}

sub add_error {
    my ( $self, $ET ) = @_;
    no warnings 'once';
    bless [
        (
            map {
                log sum pairwise { defined($a) ? $a * exp($b) : () } @$_, @$self
            } @$ET
        ),
        0
      ],
      ref($self);
}

1;
