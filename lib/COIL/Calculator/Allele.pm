package COIL::Calculator::Allele;

use strict;
use warnings;

use List::Util qw/ sum /;
use List::MoreUtils qw/ pairwise /;

use Params::Validate;
use COIL::Validate qw/ :val :grab /;

use parent 'COIL::Calculator';

=head1 NAME

COIL::Calculator::Allele - allelic likelihoods

=head1 SYNOPSIS

Compute log likelihoods for a barcode assuming independent alleles.

    my $calc = COIL::Calculator::Allele->new_from_tally( $CTA );
    my $calc_err = $calc->add_error( $e );
    my $likelihoods = $calc_err->numerics_likelihoods( $numerics );

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

	my $calc = COIL::Calculator::Allele->new_from_tally( $tally );
	my $calc = COIL::Calculator::Allele->new_from_tally(
	   $tally,
	   {
	       max_COI => 5,
	       padding => 1
	   }
	);

Create the likelihood calculation object. The unit of which is the log of the
following vector:

        / P(A|C) \
    L = | P(a|C) |
        \ P(N|C) /

Where A, a and N are the events that we observe the major allele, minor allele
or a het call at a particular assay, respectively. There is a fourth element,
log(P(X|C))=0, where X is the event that the assay failed, since failed assays
do not tell us any information about COI. To compute these elements, let R be a
count of ref calls in the different strains.

    R ~ Binomial(C, p), p = allele frequency of ref

    A = I(R=C)
    a = I(R=0)
    N = I(R>0, R<C)

    P(A) = P(R=C) = p^C
    P(a) = P(R=0) = (1-p)^C
    P(N) = 1 - P(A) - P(a)

Since we are doing everything on a log scale, the first two equations become:

    log P(A) = C * log(p)
    log P(a) = C * log(1-p)

=head2 add_error

    my $calc_err = $calc->add_error();
    my $calc_err = $calc->add_error( $e ); # default e=.05
    my $calc_err = $calc->add_error( [ $e0, $e+, $e- ])

One way to view how this is done is to do a matrix multiplication, in which the
error rates specify a stochastic matrix. Given the previous calculator vector,
we want to find

         / P(A*|C) \
    L* = | P(a*|C) |
         \ P(N*|C) /

where A*, a* and N* are the observed observations for assay. The error rate is
a stochastic matrix

    E[i][j] = P(G_i*=j|G_i=i)

        / P(A*|A), P(a*|A), P(N*|A) \
    E = | P(A*|a), P(a*|a), P(N*|a) |
        \ P(A*|N), P(a*|N), P(N*|N) /

        / 1-e, e/2, e/2 \
      = | e/2, 1-e, e/2 |, simple error rate case
        \ e/2, e/2, 1-e /

        / 1-e0-e+, e0,      e+   \
      = | e0     , 1-e0-e+, e+   |, a bit more complicated
        \ e-/2   , e-/2,    1-e- /

then the new calculator vector is the product:

    L*' = L'E = E'L

=cut

=head1 METHODS

=cut

=head2 random_numeric

    my $numeric = $calc->random_numeric( $COI );

Generate a random numeric at a particular COI.

=cut

sub random_numeric {
    my $self = shift;
    my ($COI) = validate_pos( @_, $VAL_POS_INT );
    $self->[ $COI - 1 ]->random_numeric();
}

=head2 write

    $calc->write( $fh, $digits);

=cut

sub write {
    my $self = shift;
    my $fh = grab_fh( \@_, '>' );
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

package COIL::Calculator::Allele::Level;

# The overall likelihood object stores the likelihoods at several COIs. Each
# one is a further object whose main purpose is to propagate and collate calls
# to each individual likelihood unit (below).

use strict;
use warnings;

our @ISA = ('COIL::Calculator::Level');

use List::Util 'sum';
use List::MoreUtils 'pairwise';

sub _new_from_ps {
    my $class = shift;
    ( my $unit_class = $class ) =~ s/Level$/Unit/;

    bless [ map { $unit_class->_new_from_p($_) } @{ $_[0] } ], $class;
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

package COIL::Calculator::Allele::Unit;

use strict;
use warnings;

our @ISA = ('COIL::Calculator::Unit');

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

use overload '""' => \&to_string;

sub to_string {
    join( ':', map { exp($_) } @{ $_[0] }[ 0 .. 2 ] );
}

1;
