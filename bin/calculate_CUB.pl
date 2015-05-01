#!/usr/bin/env perl

use strict;
use warnings;
use lib '../lib';
use Bio::CUA::CUB::Calculator;
use Bio::CUA::CodonTable;
use Getopt::Long;

my @args = @ARGV;
my $sep = "\t";
my $seqIO_pkg;

BEGIN{
	eval { require Bio::SeqIO; };
	
	if($@)
	{
		require Bio::CUA::SeqIO;
		$seqIO_pkg = 'Bio::CUA::SeqIO';
	}else
	{
		$seqIO_pkg = 'Bio::SeqIO';
	}
}

my $seqFile;
my $gcId;
my $outFile;
my $taiFile;
my $caiFile;
my $encMethods;
my $baseCompFile;

GetOptions(
	'gc-id:i'     => \$gcId,
	'tai-param:s' => \$taiFile,
	'cai-param:s' => \$caiFile,
	'enc-methods:s' => \$encMethods,
	'out:s' => \$outFile,
	'base-comp:s'	=> \$baseCompFile
);

$seqFile = shift or &usage();

$gcId ||= 1;
$outFile ||= '-';
$encMethods ||= 'enc';

my @encs = split ',', $encMethods;

my %baseCompositions;
if(grep /^encp/,@encs) # there are methods needing background data
{
	die "option --base-comp is missing for encp and its related methods:$!" 
	unless($baseCompFile);
	open(B,"< $baseCompFile") or die "Can not open $baseCompFile:$!";
	while(<B>)
	{
		next if /^#/;
		chomp;
		my @fields = split "\t";
		# we may implement checking the validity of the data here in
		# future
		$baseCompositions{$fields[0]} = [@fields[1..4]];
	}
	close B;
}

warn "Step 1: build analyzer\n";

my $table = Bio::CUA::CodonTable->new(-id => $gcId) or die $!;
my %params; # to CUB analyzer
%params = (
	-codon_table  => $table,
);
$params{'-tai_values'} = $taiFile if($taiFile);
$params{'-cai_values'} = $caiFile if($caiFile);

my $cub  = Bio::CUA::CUB::Calculator->new(%params) or die $!;

warn "Step 2: calculating parameters\n";
my $io = $seqIO_pkg->new(-file => $seqFile, -format => 'fasta') or die
"Can't create sequence IO on $seqFile:$!";

open(O, "> $outFile") or die "Can not open $outFile:$!";
my @header = qw/seq_id length AAs GC GC3/;
push @header, @encs;
push @header, 'tai' if($taiFile);
push @header, 'cai' if($caiFile);

my @allAAtypes = sort($table->all_amino_acids); # make the amino acids
# always in alphabet order
print O "# produced by $0 @args\n";
print O "# amino acid order for AAs: ",join(',',@allAAtypes),"\n";
print O '#', join($sep, @header),"\n";
my $counter = 0;
while(my $seq = $io->next_seq)
{
	my $codonList = $cub->get_codon_list($seq) or 
	(warn "Get codons for ".$seq->id." failed\n" and next);
	my %AAs;
	my $totalPepLen = 0;
	while(my ($codon, $cnt) = each %$codonList)
	{
		next if($table->is_stop_codon($codon));
		my $AA = $table->translate($codon) or 
		die "Can not translate $codon:$!";
		$AAs{$AA} += $cnt;
		$totalPepLen += $cnt;
	}

	my $baseComp = $baseCompositions{$seq->id} || 'NA'
	if(%baseCompositions);
	# call (codons, minTotal, base composition)
	# set undef value to NA
	my @encResult = map { $cub->$_($codonList,undef,$baseComp) || 'NA' } @encs;

	my $aaCounts = join(',', map { $AAs{$_} || 0 } @allAAtypes);
	my $gcFrac = $cub->gc_frac($codonList);
	my $gc3Frac = $cub->gc_frac($codonList,3);
	my $tai = $cub->tai($codonList) or die $! if($taiFile);
	my $cai = $cub->cai($codonList) or die $! if($caiFile);
	my @result = ($seq->id, $totalPepLen, $aaCounts);
	# set ENC which is > 62 to 62
	push @result, map { $_ ne 'NA' and $_ > 62? 62 : _format_num($_) } 
	     ($gcFrac, $gc3Frac, @encResult);
	push @result, _format_num($tai) if($tai);
	push @result, _format_num($cai) if($cai);

	print O join($sep, @result),"\n";

	warn "# $counter sequences have been analyzed\n" 
	if(++$counter %	500 == 0);
}

close O;

warn "Work [$counter sequences] is done!!\n";

exit 0;

sub _format_num
{
	my ($num) = @_;
	return $num if($num eq 'NA');
	my $field = 7;
	return sprintf("%.*f",, $field, $num);
}

sub usage
{
	print <<USAGE;
Usage: $0 [options] <seq-file>

This program can read into sequences from fasta-formated file
<seq-file> and output information of each sequence, including amino
acid composition, codon usage bias (ENC, CAI, tAI, Fop, etc), GC
content.

Options:

--gc-id:  genetic code table ID used for deriving amino acids for
sequences in input file. Default is 1, i.e., standard table

--tai-param:  the file containing tAI values for each codon in the
format 'codon<tab>tAI_value'. If not provided, tAI values would not be
computed.

--cai-param: similar to --tai-param, except for CAI values.

--enc-methods: what methods to be used to calculate ENC (codon usage
bias). Available values are enc, enc_r, encp, and encp_r, Check module
Bio::CUA::CUB::Calculator to see what each method does. Default is
enc. One can specify more than one method such as 'enc,encp,enc_r'.

--base-comp: background base compositions in the format of 
seq_id	#A	#T	#C	#G
where #A/#T/#C/#G are counts or fractions of bases in background data
such as introns. This option is mandatory if --enc-methods has
arguments 'encp' or 'encp_r'. One line per sequence. For sequences
mising this information, an 'NA' will be returned for corresponding enc
values.

--out:  the file to store the results. Default is to standard output.

USAGE
	
	exit 1;
}
