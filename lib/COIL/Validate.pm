package COIL::Validate;

use strict;
use warnings;
no warnings 'once';

use Params::Validate;
use Carp;
use Scalar::Util qw/ looks_like_number reftype /;
use List::Util qw/ sum /;

# (1) INITIALIZE EXPORTING
use Exporter 'import';

# regular expressions
our @EXPORT_RE = qw/
  $RE_NUC
  $RE_SALLELE
  $RE_STRAIN

  $RE_UNIT
  $RE_ALLELE
  $RE_BARCODE

  $RE_NUM
  $RE_NALLELE
  $RE_NUMERIC
  /;

# validations for use in Params::Validate
our @EXPORT_VAL = qw/
  $VAL_NUM
  $VAL_POS_INT
  $VAL_NON_NEG_INT
  $VAL_POS_REAL
  $VAL_NON_NEG_REAL
  $VAL_NON_POS_REAL
  $VAL_PROB
  $VAL_PADDING

  $VAL_NUMS
  $VAL_PROBS
  $VAL_NON_POS_REALS
  $VAL_DIST
  $VAL_LOG_DIST
  $VAL_DISTS
  $VAL_LOG_DISTS

  $VAL_NUC
  $VAL_ALLELE
  $VAL_NALLELE
  $VAL_SALLELE

  $VAL_BARCODE_STR
  $VAL_STRAIN_STR
  $VAL_NUMERIC_STR

  $VAL_BARCODE
  $VAL_STRAIN
  $VAL_NUMERIC

  $VAL_BARCODES
  $VAL_STRAINS
  $VAL_NUMERICS
  /;

# functions to do standalone validation
our @EXPORT_FUN = qw/
  val_allele
  val_sallele
  val_nallele

  val_barcode
  val_strain
  val_numeric

  val_barcodes
  val_strains
  val_numerics

  val_arrayrefs_same_length
  /;

# these grab something from @_ before validation
our @EXPORT_GRAB = qw/
  grab_fh
  /;

our @EXPORT_OK = ( @EXPORT_RE, @EXPORT_VAL, @EXPORT_FUN, @EXPORT_GRAB );

our %EXPORT_TAGS = (
    all  => \@EXPORT_OK,
    re   => \@EXPORT_RE,
    val  => \@EXPORT_VAL,
    fun  => \@EXPORT_FUN,
    grab => \@EXPORT_GRAB
);

# (2) GENERATORS
sub _gen_VAL_NUM {
    return {
        type     => Params::Validate::SCALAR,
        callback => {
            'looks like a number' => sub { looks_like_number( $_[0] ) },
            @_
        }
    };
}

# the following functions work in tandem
# _gen_val_mult creates a function that evaluates every entry within an array
# _gen_val_arrayref checks that its an arrayref, then calls the mult function
# _gen_VAL_ARRAYREF does the same but for a Params::Validate

sub _gen_val_mult {
    my ($val_fun) = @_;

    return sub {
        foreach my $e ( @{ $_[0] } ) { return unless &$val_fun($e) }
        return 1;
    };
}

sub _gen_val_arrayref {
    my ($val_fun) = @_;

    return sub {
        return unless ( reftype( $_[0] ) eq 'ARRAY' );
        return unless &$val_fun( $_[0] );
        return 1;
    };
}

sub _gen_val_arrayrefs {
    my ($val_fun) = @_;

    return sub {
        return unless ( reftype( $_[0] ) eq 'ARRAY' );
        return unless &$val_fun( $_[0] );
        return unless val_arrayrefs_same_length( $_[0] );
        return 1;
    };
}

sub _gen_VAL_ARRAYREF {
    my ( $name, $val_fun ) = @_;

    return {
        type      => Params::Validate::ARRAYREF,
        callbacks => { "contains valid $name" => $val_fun }
    };
}

sub _gen_VAL_ARRAYREFS {
    my ( $name, $val_fun ) = @_;

    return {
        type      => Params::Validate::ARRAYREF,
        callbacks => {
            "contains valid $name" => $val_fun,
            "$name same length"    => \&val_arrayrefs_same_length
        }
    };
}

sub val_arrayrefs_same_length {
    my $l = scalar( @{ $_[0][0] } );
    for ( my $i = 1 ; $i < @{ $_[0] } ; $i++ ) {
        return if ( $l != scalar( @{ $_[0][1] } ) );
    }
    return 1;
}

sub approx {
    my ( $got, $exp, $margin ) = @_;
    $margin ||= .01;
    ( $got >= ( 1 - $margin ) * $exp ) && ( $got <= ( 1 + $margin ) * $exp );
}

# (3) VALIDATIONS FOR NUMBERS
our $VAL_NUM          = _gen_VAL_NUM();
our $VAL_POS_INT      = { regex => qr/^[1-9]\d*$/ };
our $VAL_NON_NEG_INT  = { regex => qr/^\d+$/ };
our $VAL_POS_REAL     = _gen_VAL_NUM( '>0' => sub { ( $_[0] > 0 ) } );
our $VAL_NON_NEG_REAL = _gen_VAL_NUM( '>=0' => sub { ( $_[0] >= 0 ) } );
our $VAL_NON_POS_REAL = _gen_VAL_NUM( '<=0' => sub { ( $_[0] <= 0 ) } );
our $VAL_PROB         = _gen_VAL_NUM(
    '>=0' => sub { ( $_[0] >= 0 ) },
    '<=1' => sub { ( $_[0] <= 1 ) }
);

our $VAL_PADDING = { %$VAL_NON_NEG_REAL, default => 1 };

*_val_nums = _gen_val_mult( \&looks_like_number );
*val_nums  = _gen_val_arrayref( \&_val_nums );
our $VAL_NUMS = _gen_VAL_ARRAYREF( \&_val_nums );

sub val_prob { looks_like_number( $_[0] ) && ( $_[0] >= 0 ) && ( $_[0] <= 1 ) }
sub val_non_pos_real { looks_like_number( $_[0] ) && ( $_[0] <= 0 ) }

*_val_probs         = _gen_val_mult( \&val_prob );
*_val_non_pos_reals = _gen_val_mult( \&val_non_pos_real );

sub val_sum_1 { approx( sum( @{ $_[0] } ), 1 ) }

sub val_exp_sum_1 {
    approx( sum( map { exp } @{ $_[0] } ), 1 );
}

*val_probs         = _gen_val_arrayref( \&_val_probs );
*val_non_pos_reals = _gen_val_arrayref( \&_val_non_pos_reals );

sub val_dist     { val_probs(@_)         && val_sum_1(@_) }
sub val_log_dist { val_non_pos_reals(@_) && val_exp_sum_1(@_) }

*_val_dists     = _gen_val_mult( \&val_dist );
*_val_log_dists = _gen_val_mult( \&val_log_dist );

our $VAL_PROBS = _gen_VAL_ARRAYREF( 'probabilities', \&_val_probs );
our $VAL_NON_POS_REALS =
  _gen_VAL_ARRAYREF( 'non positive reals' => \&_val_non_pos_reals );

our $VAL_DIST = _gen_VAL_ARRAYREF( 'probabilities', \&_val_probs );
$VAL_DIST->{callbacks}{'sum to 1'} = \&val_sum_1;

our $VAL_LOG_DIST =
  _gen_VAL_ARRAYREF( 'log probabilities', \&_val_non_pos_reals );
$VAL_LOG_DIST->{callbacks}{'sum to 1'} = \&val_exp_sum_1;

our $VAL_DISTS = _gen_VAL_ARRAYREFS( distributions => \&_val_dists );
our $VAL_LOG_DISTS =
  _gen_VAL_ARRAYREFS( 'log distributions' => \&_val_log_dists );

# (4) STRING VALIDATION
our $RE_NUC     = qr/[ACGT]/;
our $RE_SALLELE = qr/^$RE_NUC$/;
our $RE_STRAIN  = qr/^$RE_NUC+$/;

our $RE_UNIT    = qr/[ACGTNX]/;
our $RE_ALLELE  = qr/^$RE_UNIT$/;
our $RE_BARCODE = qr/^$RE_UNIT+$/;

our $RE_NUM     = qr/[0-3]/;
our $RE_NALLELE = qr/^$RE_NUM$/;
our $RE_NUMERIC = qr/^$RE_NUM+$/;

our $VAL_ALLELE  = { regex => $RE_ALLELE };
our $VAL_SALLELE = { regex => $RE_SALLELE };
our $VAL_NALLELE = { regex => $RE_NALLELE };

our $VAL_BARCODE_STR = { regex => $RE_BARCODE };
our $VAL_STRAIN_STR  = { regex => $RE_STRAIN };
our $VAL_NUMERIC_STR = { regex => $RE_NUMERIC };

# validate barcode-type object
our $VAL_BARCODE = _gen_VAL_ARRAYREF( alleles => \&_val_alleles );
our $VAL_STRAIN  = _gen_VAL_ARRAYREF( alleles => \&_val_salleles );
our $VAL_NUMERIC = _gen_VAL_ARRAYREF( alleles => \&_val_nalleles );

# validate arrays of barcode-type objects
our $VAL_BARCODES = _gen_VAL_ARRAYREFS( 'barcodes', \&_val_barcodes );
our $VAL_STRAINS  = _gen_VAL_ARRAYREFS( 'strains',  \&_val_strains );
our $VAL_NUMERICS = _gen_VAL_ARRAYREFS( 'numerics', \&_val_numerics );

# functions for validating individual alleles
sub val_allele  { $_[0] =~ $RE_ALLELE }
sub val_sallele { $_[0] =~ $RE_SALLELE }
sub val_nallele { $_[0] =~ $RE_NALLELE }

# and their plurals
*_val_alleles  = _gen_val_mult( \&val_allele );
*_val_salleles = _gen_val_mult( \&val_sallele );
*_val_nalleles = _gen_val_mult( \&val_nallele );

# functions for validating objects
*val_barcode = _gen_val_arrayref( \&_val_alleles );
*val_strain  = _gen_val_arrayref( \&_val_salleles );
*val_numeric = _gen_val_arrayref( \&_val_nalleles );

# and their plurals
*_val_barcodes = _gen_val_mult( \&val_barcode );
*_val_strains  = _gen_val_mult( \&val_strain );
*_val_numerics = _gen_val_mult( \&val_numeric );

# functions to validate arrays of barcode-type objects
*val_barcodes = _gen_val_arrayrefs( \&_val_barcodes );
*val_strains  = _gen_val_arrayrefs( \&_val_strains );
*val_numerics = _gen_val_arrayrefs( \&_val_numerics );

# GRAB FUNCTIONS

# helper function that opens file handle and spits out error messages

sub grab_fh {
    my ( $params, $mode ) = @_;
    my $long_mode;

    $mode ||= '<';
    if    ( $mode eq '>' ) { $long_mode = 'writing' }
    elsif ( $mode eq '<' ) { $long_mode = 'reading' }
    else                   { croak qq{Invalid mode "$mode"} }

    my $file = shift(@$params);
    unless ($file) { return $mode eq '<' ? *STDIN : *STDOUT }

    my $type = ref($file) || ref( \$file );

    if ( $type eq 'SCALAR' ) {
        open my $fh, ( $mode || '<' ), $file
          or croak qq{Unable to open file "$file" for $long_mode};
        return $fh;
    }
    if ( $type eq 'GLOB' ) { return $file }

    croak qq{Invalid type "$type"};
}

1;
