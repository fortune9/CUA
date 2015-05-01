#!/usr/bin/env perl

=pod

a pod test in script

=cut

use strict;
use warnings;
use lib '../lib';
use Bio::CUA::CUB::Builder;
use Bio::CUA::CodonTable;
use Getopt::Long;
use Fcntl qw/:seek/;
use File::Sort qw/sort_file/;
use File::Temp qw/ tempfile tempdir /;

my $seqIO_pkg;

BEGIN{
	eval { require Bio::SeqIO; };

	if($@) # bioperl is not installed
	{
		require Bio::CUA::SeqIO;
		$seqIO_pkg = 'Bio::CUA::SeqIO';
	}else
	{
		$seqIO_pkg = 'Bio::SeqIO';
	}
}

my $sep = "\t"; # field separator
my $seqFile;
my $expFile;
my $gcId;
my $outFile;
my $select;
my $background;
my $normMethod;
my $minTotal; # the minimal count of an amino acid for being
# considered to calculate CAI parameters

GetOptions(
	'seq-file=s'  => \$seqFile,
	'exp-file:s'  => \$expFile,
	'select:s'  => \$select,
	'background:s'  => \$background,
	'out-file:s'  => \$outFile,
	'gc-id:i'  => \$gcId,
	'method:s' => \$normMethod
);

&usage() unless($seqFile);

if($select)
{
	die "option '--exp-file' is needed as option '--select' is provided:$!" 
	unless($expFile);
}

$gcId ||= 1;
$outFile ||= '-';
$select ||= 'all';
$select = lc($select);
$normMethod ||= 'max';

warn "# Step 1: tabulate codons from sequences\n";
my $table = Bio::CUA::CodonTable->new(-id => $gcId) 
   or die "Creating codon table failed:$!";
my $builder = Bio::CUA::CUB::Builder->new(
	           -codon_table  => $table
		   )
   or die "Creating analyzer failed:$!";

my $codonListRef;
my ($sortedFh, $sortedExpFile);
if($expFile)
{
	# sort this $expFile first
	($sortedFh, $sortedExpFile) = tempfile();
	#warn "$sortedExpFile\n";
	sort_file({
		I => $expFile,
		o => $sortedExpFile,
		k => '2,2n',
		t => "\t" # field separator
	});
	my $idsRef;
	if($select eq 'all')
	{
		$idsRef = choose_ids($sortedFh, 'all');
	}elsif($select < 1) # a fraction
	{
		my $total = num_of_lines($expFile);
		$select = int($total * $select);
		# largest is at the end, -r option for sort does not work
		$idsRef = choose_ids($sortedFh, 'high', $select);
	}elsif($select > 1) # number of genes
	{
		$select = int($select);
		$idsRef = choose_ids($sortedFh, 'high', $select);
	}else
	{
		die "Unknown option '$select' for --select:$!";
	}
	
	$codonListRef = &read_codons($seqFile,$idsRef);
	# output_ids($idsRef,'positive');
}else # consider codons in all input sequences
{
	$codonListRef = &read_codons($seqFile);
}

# now consider background data

my $backCodonListRef;
if($background)
{
	if(is_numeric($background))
	{
		die "option '--exp-file' is needed for selecting background data:$!" 
		unless($sortedFh);

		my $idsRef;
		if($background < 1) # a fraction
		{
			my $total = num_of_lines($expFile);
			$background = int($total * $background);
			# largest is at the end, -r option for sort does not work
			$idsRef = choose_ids($sortedFh, 'low', $background);
		}elsif($background > 1) # number of genes
		{
			$background = int($background);
			$idsRef = choose_ids($sortedFh, 'low', $background);
		}else
		{
			die "Unknown option '$background' for --background:$!";
		}
	
		$backCodonListRef = &read_codons($seqFile,$idsRef);

		#output_ids($idsRef,'negative');
	}elsif(-f $background) # a sequence file
	{
		$backCodonListRef = &read_codons($background) or 
		die "reading codons from $background failed:$!";
	}else
	{
		die "Unknown data '$background' for option --background:$!";
	}
}

warn "# Step 2: calculate CAIs for codons\n";
if($background)	
{
	$builder->build_b_cai($codonListRef, $backCodonListRef,
		$minTotal, $outFile) or die "building codons' b_CAI failed:$!";
}else
{
	$builder->build_cai($codonListRef, $minTotal, $normMethod,
		$outFile) or 
	die "building codons' CAI failed:$!";
}

# remove temporary files
if($sortedExpFile)
{
	close $sortedFh;
	unlink $sortedExpFile;
}

warn "Work done!!\n";

exit 0;

sub output_ids
{
	my ($hash, $name) = @_;

	open(O,"> tmp_$name") or die $!;
	print O join("\n", keys %$hash), "\n";
	close O;
}

sub is_numeric
{
	my $num = shift;
	if($num =~ /^[+-]?[eE\d\.\-]+$/)
	{
		return 1;
	}else
	{
		return 0;
	}
}

# get the number of lines in a file
sub num_of_lines
{
	my ($file) = @_;
	my $totalNum;
	if(ref($file) eq 'GLOB')
	{
		$. = 0; # reset number to 0
		my $currPos = tell($file);
		seek($file, 0, SEEK_SET); # set to beginning
		1 while(<$file>);
		$totalNum = $.; # a special variable
		seek($file,$currPos,SEEK_SET); # set back the pointer
	}else
	{
		my $fh;
		open($fh, "< $file") or die "Can not open $file:$!";
		1 while(<$fh>);
		$totalNum = $.; # a special variable
		close $fh;
	}
	return $totalNum;
}

# read codons from given sequence IDs
sub read_codons
{
	my ($seqFile,$idsRef) = @_;
	# if $idsRef is not a hash reference, all the sequences in the
	# input file will be scanned

	if(!defined($idsRef)) # no selected IDs
	{
		return $builder->get_codon_list($seqFile);
	}

	my $io = $seqIO_pkg->new(-file => $seqFile, -format => 'fasta');
	my %codonList; # store all the codons of selected sequences
	my $seqCnt = 0;

	while(my $seq = $io->next_seq)
	{
		next unless($idsRef->{$seq->id});
		my $localCodons = $builder->get_codon_list($seq) or 
		(warn "Counting codons in ".$seq->id." failed\n" and next);
		while(my ($codon, $cnt) = each %$localCodons)
		{
			$codonList{$codon} += $cnt;
		}
		$seqCnt++;
	}

	return \%codonList;
}

# choose IDs from an input file
sub choose_ids
{
	my ($fh, $method, $cnt) = @_;
	$method = lc($method);
	# set the filehandle to the beginning
	seek($fh,0,SEEK_SET);
	my %data;
	my $counter = 0;
	if($method eq 'all')
	{
		while(<$fh>)
		{
			chomp;
			my @fields = split $sep;
			$data{$fields[0]}++;
		}
		# return _read_column($fh, 1);
	}elsif($method eq 'low')
	{
		while(<$fh>)
		{
			chomp;
			my @fields = split $sep;
			$data{$fields[0]}++;
			last unless(++$counter < $cnt);
		}
	}elsif($method eq 'high')
	{
		my $toSkip = num_of_lines($fh) - $cnt;
		while(<$fh>)
		{
			next unless(++$counter > $toSkip);
			chomp;
			my @fields = split $sep;
			$data{$fields[0]}++;
		}
	}else
	{
		die "Unknown method $method for choosing IDs in 'choose_ids':$!";
	}

	return \%data;
}

sub _read_column
{
	# the $fh may not be at the beginning
	my ($fh, $colNum, $cnt) = @_;
	my %data;
	while(<$fh>) 
	{
		last if(defined($cnt) and $cnt == 0);
		chomp;
		my @fields = split $sep;
		$data{$fields[$colNum - 1]}++; 
	}

	return \%data;
}

sub usage
{
	print <<USAGE;
Usage: $0 [options]

This program reads into fasta-formatted sequences and ouput CAI values
for each codon. These input sequences are supposed to be highly
expressed genes in an organism and thus presumably use efficient codons.

Mandatory options:

--seq-file:  a file containing sequences in fasta format

Auxiliary options:

--exp-file:  a file containing a list of sequence IDs with one ID per
line, and optionally each ID followed by sequence expression levels.
The two fields are separated by <tab>. See option --select how this
file affects the program's process.

--select:  determine how sequences are selected from sequence file.
Available values are as follows:
	all  = all sequences with IDs in the expression file given by --exp-file
	0.## = a fraction such as 0.30, then top 30% sequences in the
	expression file ranked after gene expression (from high to low) are 
	selected for CAI calculation.
	
	###  = an integer such as 200, then the top 200 sequences (based on
	gene expression rank) in the expression file are selected for CAI
	calculation
	
	Default is 'all'.
	Note: the sequence IDs must match sequence IDs in the sequence file to
	be able to select those sequences.

--background: specify the background data such as lowly expressed
genes from which the background codon usage is derived. This argument
is optional, but if provided, 'background-normalization' method is
used, see below for details.
Acceptable values are as follows:
	0.## = a fraction such as 0.30, then bottom 30% sequences in the
	expression file ranked after gene expression (from high to low).
	
	###  = an integer such as 200, then the bottom 200 sequences (based on
	#gene expression rank) in the expression file are selected.
	filename = a file containing protein-coding sequences which will be
	parsed for background codon usage.

--gc-id:  id of genetic code table. See NCBI genetic code
http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes
for valid IDs. Default is 1, i.e., standard genetic code.

--method: the method to normalize CAI values of synonymous codons of
the same amino acid. Default is 'max'; the alternatives are 'mean' or
'background-normalization'. 
The 'max' method is the one used by <Sharp and Li, 1987, NAR>. 
'mean' divides each codon's RSCU by its expectation under even usage 
of synonymous codons. For example, for an amino acid with four synymous 
codons, the expectation is 0.25 for each synonymous codon.
'background-normalization' divides each codon's RSCU by the
corresponding RSCU from backgound sequence data, then these normalized
RSCU values follows the 'max' method to derive CAI values.

--out-file: the file to store the result. Default is standard output.

Author:  Zhenguo Zhang
Contact: zhangz.sci\@gmail.com
Created:
Wed Mar  4 16:42:24 EST 2015

USAGE
	
	exit 1;
}
