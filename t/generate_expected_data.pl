#!/usr/bin/env perl

use strict;
use lib '/home/zzhang65/Scripts/Perl/CPAN_dist/Bio-CUA/lib';
use Bio::CUA::CodonTable;
use Bio::CUA::Summarizer;
use Bio::CUA::SeqIO;
use Bio::CUA::CUB::Calculator;
use Bio::CUA::CUB::Builder;

warn "The program is to generate results for test purpose\n";
warn "Shall we continue?[Y/N]\n";

my $in = <STDIN>;
unless($in =~ /^Y/i)
{
	warn "Stopped\n";
	exit 0;
}

my $seqFile = 'test.fa'; # sequence data
my $tRNAFile = 'dmel_r5_apr2006.tRNA_copy'; # tRNA copy number

#my $outFile = shift || '-';

my $builder = Bio::CUA::CUB::Builder->new( -codon_table => 1 );

warn "# Generating expected results for Bio-CUA\n";

# 1. codon table part
my $codonTable = $builder->codon_table;
my $sect = 'Bio::CUA::CodonTable';
print ">>$sect\n";
print ">>>codon_to_AA\n";
my $codon2AA = $codonTable->codon_to_AA_map or die $!;
foreach my $c (sort keys(%$codon2AA))
{
	printf("%s\t%s\n", $c, $codon2AA->{$c});
}
print ">>>degeneracy\n";
my $codonDegeneracy = $codonTable->codon_degeneracy or die $!;
foreach my $deg (sort keys %$codonDegeneracy)
{
	foreach my $aa (keys %{$codonDegeneracy->{$deg}})
	{
		my $codons = $codonDegeneracy->{$deg}->{$aa};
		printf("%d\t%s\t%s\n", $deg, $aa, join(',', @$codons));
	}
}
print "<<$sect\n\n";


# 2. SeqIO part
$sect = 'Bio::CUA::SeqIO';
my $io = $sect->new(-file => $seqFile);
print ">>$sect\n";
print ">>>length\n";
while(my $seq = $io->next_seq)
{
	printf("%s\t%d\n", $seq->id, $seq->length);
}
print "<<$sect\n\n";

# 3. summarizer part
$sect = 'Bio::CUA::Summarizer';
my $sum = $sect->new(-codon_table => 1);
my $codonList = $sum->tabulate_codons($seqFile) or die $!;
print ">>$sect\n";
print ">>>codon_count\n";
_write_hash($codonList);
print "<<$sect\n\n";

# 4. builder part
$sect = 'Bio::CUA::CUB::Builder';
# traditional method, average, and background
my $rscu = $builder->build_rscu($seqFile,5,0.5); 
my $caiMax  = $builder->build_cai($seqFile,5,'max'); 
my $caiMean = $builder->build_cai($seqFile,5,'mean');
my $caiBack = $builder->build_b_cai($seqFile,$seqFile,5);
my $tai = $builder->build_tai($tRNAFile);
print ">>$sect\n";
print ">>>rscu\n";
_write_hash($rscu);
print ">>>cai_max\n";
_write_hash($caiMax);
print ">>>cai_mean\n";
_write_hash($caiMean);
print ">>>cai_back\n";
_write_hash($caiBack);
print ">>>tai\n";
_write_hash($tai);
print "<<$sect\n\n";

# 5. calculator part
$sect = 'Bio::CUA::CUB::Calculator';
my $calc = $sect->new(
	          -codon_table => 1,
			  -CAI_values  => $caiMax,
			  -tAI_values  => $tai
		  );

# calculate CAI/tAI/ENC/ENC_r for each sequence
$io = Bio::CUA::SeqIO->new(-file => $seqFile);
my %seqCAIs;
my %seqtAIs;
my %seqENCs;
my %seqENC_rs;

while(my $seq = $io->next_seq)
{
	$seqCAIs{$seq->id} = $calc->cai($seq);
	$seqtAIs{$seq->id} = $calc->tai($seq);
	$seqENCs{$seq->id} = $calc->enc($seq,5);
	$seqENC_rs{$seq->id} = $calc->enc_r($seq,5);
}

print ">>$sect\n";
print ">>>seq_cai\n";
_write_hash(\%seqCAIs);
print ">>>seq_tai\n";
_write_hash(\%seqtAIs);
print ">>>seq_enc\n";
_write_hash(\%seqENCs);
print ">>>seq_enc_r\n";
_write_hash(\%seqENC_rs);
print "<<$sect\n\n";


warn "Work done!!\n";

exit 0;

sub _write_hash
{
	my $hashRef = shift;
	foreach (sort keys(%$hashRef))
	{
		printf("%s\t%s\n",$_, $hashRef->{$_});
	}

	return 1;
}

