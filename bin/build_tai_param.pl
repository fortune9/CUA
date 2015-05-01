#!/usr/bin/env perl

use strict;
use warnings;
use lib '../lib';
use Bio::CUA::CUB::Builder;
use Bio::CUA::CodonTable;

my $tRNAFile = shift;
my $gcId = shift;
my $outFile = shift;

&usage() unless($tRNAFile);

$gcId ||= 1;
$outFile ||= '-';

my $table = Bio::CUA::CodonTable->new(-id => $gcId);
my $builder = Bio::CUA::CUB::Builder->new(
	           -codon_table  => $table
		   );

$builder->build_tai($tRNAFile, $outFile) or 
die "building tAI failed:$!";

warn "Work done!!\n";

exit 0;

sub usage
{
	print <<USAGE;
Usage: $0 <tRNA-copy-number> [<genetic-code-id> <outfile>]

USAGE
	
	exit 1;
}
