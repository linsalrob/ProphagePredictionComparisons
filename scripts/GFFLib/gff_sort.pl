#!/bin/env perl

# Sort GFF file based on the numerical part of contig/scaffolds names
# and then by the start coord. of genes

use strict;

use FindBin;
use lib "$FindBin::Bin";

use GFFFile;


my $usage =
"gff_sort.pl <GFF file> <output: sorted GFF file>\n\n";

die $usage if scalar(@ARGV) != 2;

my $in_gff        = $ARGV[0];
my $out_gff       = $ARGV[1];

my $gffFile = GFFFile::new($in_gff);
$gffFile->read();

my $gffGenes = $gffFile->get_genes_hash();

# Ordering genes based on template name and start coord
my @gffGenesArray = values %{$gffGenes};

# Sort by the numerical part of chrom and then start coord of gene
@gffGenesArray =
  sort { 
  		my ($a_num_chrom) = ( $a->get_chrom() =~ /(\d+)/ );
  		my ($b_num_chrom) = ( $b->get_chrom() =~ /(\d+)/ );
  		
  		 return $a_num_chrom <=> $b_num_chrom || $a->get_start() <=> $b->get_start() }
  @gffGenesArray;


open OUT, ">$out_gff" or die "Unable to write on file $out_gff\n";
for my $currGene (@gffGenesArray) {	
	print OUT $currGene->toGFF();	
}
close(OUT);





