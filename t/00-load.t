#!perl -T
use 5.006;
use strict;
use warnings;
use Test::More;

plan tests => 3;

BEGIN {
    use_ok( 'Bio::CUA::CUB::Builder' ) || print "Bail out!\n";
    use_ok( 'Bio::CUA::CUB::Calculator' ) || print "Bail out!\n";
    use_ok( 'Bio::CUA' ) || print "Bail out!\n";
}

diag( "Testing Bio::CUA::CUB::Builder $Bio::CUA::CUB::Builder::VERSION, Perl $], $^X" );
