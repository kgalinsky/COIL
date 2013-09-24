package COIL::Likelihood::Allele;

use strict;
use warnings;

use List::Util qw/ sum /;
use List::MoreUtils qw/ pairwise /;

use Params::Validate;
use COIL::Validate ':val';

use COIL '_fh';

use base 'COIL::Likelihood';

=head1 NAME

COIL::Likelihood::Allele - allelic likelihoods

=head1 SYNOPSIS

Compute log likelihoods for a barcode assuming independent alleles.

    my $CLA = COIL::Likelihood::Allele->new_from_tally( $CTA );
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

=head2 new_from_tally

	my $CLA = COIL::Likelihood::Allele->new_from_tally( $tally );
	my $CLA = COIL::Likelihood::Allele->new_from_tally(
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

=head1 METHODS

=cut

=head2 add_error

    my $CLA_E = $CLA->add_error();
    my $CLA_E = $CLA->add_error( $e ); # default e=.05
    my $CLA_E = $CLA->add_error( [ $e0, $e+, $e- ])

See COIL::Likelihood for an explanation of the error matrix. The likelihood
with errors can be computed using matrix multiplication on the stochastic
matrix.

    L*' = L'E
    L*  = E'L = EL

=cut

=head2 random_numeric

    my $numeric = $CLA->random_numeric( $COI );

Generate a random numeric at a particular COI.

=cut

sub random_numeric {
    my $self = shift;
    my ($COI) = validate_pos( @_, $VAL_POS_INT );
    $self->[ $COI - 1 ]->random_numeric();
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

# The overall likelihood object stores the likelihoods at several COIs. Each
# one is a further object whose main purpose is to propagate and collate calls
# to each individual likelihood unit (below).

use strict;
use warnings;

use base 'COIL::Likelihood::Level';

use List::Util 'sum';
use List::MoreUtils 'pairwise';

sub _new_from_ps {
    bless [ map { COIL::Likelihood::Allele::Unit->_new_from_p($_) }
          @{ $_[1] } ], $_[0];
}

# log L(C|G) = sum log L(C|G_i)
sub _numeric_likelihood {
    my ( $self, $numeric ) = @_;
    no warnings 'once';
    sum pairwise { $a->[$b] } @$self, @$numeric;
}

# Random G is collection of random G_is
sub random_numeric {
    [ map { $_->random_numeric } @{ $_[0] } ];
}

package COIL::Likelihood::Allele::Unit;

use strict;
use warnings;

use base 'COIL::Likelihood::Unit';

use List::Util 'sum';
use List::MoreUtils 'pairwise';

use Params::Validate;
use COIL::Validate ':val';

sub _new_from_tally { $_[0]->_new_from_p( $_[1]->p( $_[2] ) ) }

# create a new object from a probability
sub new_from_p {
    shift->_new_from_p( validate_pos( @_, $VAL_PROB ) );
}

sub _new_from_p {
    bless [ log( $_[1] ), log( 1 - $_[1] ), '-inf', 0 ], $_[0];
}

# increment to the next COI
sub _increment {
    my $obj = bless [ $_[0][0] + $_[1][0], $_[0][1] + $_[1][1], undef, 0 ],
      ref( $_[0] );
    $obj->[2] =
      log( 1 - exp( $obj->[0] ) - exp( $obj->[1] ) );
    return $obj;
}

# perform the cross product
sub _add_error {
    my ( $self, $E ) = @_;
    bless [
        (
            map {
                log(
                    exp( $self->[0] ) * $E->[0][$_] +
                      exp( $self->[1] ) * $E->[1][$_] +
                      exp( $self->[2] ) * $E->[2][$_] )
            } ( 0 .. 2 )
        ),
        0
      ],
      ref($self);
}

# Random G_i
sub random_numeric {
    my $r = rand(1);

    for ( my $n = 0 ; $n < 3 ; $n++ ) {
        my $p = exp( $_[0][$n] );
        return $n if ( $r < $p );
        $r -= $p;
    }

    croak(
        sprintf(
            'Random value greater than probabilities: %g %g %g',
            exp( $_[0][0] ),
            exp( $_[0][1] ),
            exp( $_[0][2] )
        )
    );
}

1;
