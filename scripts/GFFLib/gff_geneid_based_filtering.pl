#!/bin/env perl
use strict;

use FindBin;
use lib "$FindBin::Bin";

use GFFFile;

my $usage = "gff_geneid_based_filtering.pl <id list file> <GFF file> <GFF out> [remove]\n\n";

die $usage if ( scalar( @ARGV ) != 3 && scalar( @ARGV ) != 4 );

my $idFile = $ARGV[0];
my $outFile = $ARGV[2];

# Indicate if the gene list should be either kept or removed from the GFF file 
my $remove = 0;
$remove = 1 if( defined( $ARGV[3] ) );

# Genes list
my %gene_list;

my $num_genes_list = 0;
my %found;

# Read file listing genes
open IN, $idFile or die "Unable to open file $idFile\n"; 
while( <IN> ){
	my $line = $_;
	chomp $line;
	$gene_list{ $line } = 1;
	$num_genes_list++;
	$found{ $line } = -1;
}
close(IN);

# Reading GFF file	
my $gffFile = GFFFile::new( $ARGV[1] );
$gffFile->read();

my $gffGenes = $gffFile->get_genes_hash();



open OUT, ">$outFile" or die "Unable to open file $outFile\n"; 
my $num_genes_found = 0;
# Printing only genes that should be kept
for my $currGene (values %{$gffGenes} ){
	my $gene_id = $currGene->get_id();
	if ( $gene_list{ $gene_id } == 1 ){
		$found{$gene_id} =  1;
		$num_genes_found++;	
	}
	if( $remove ){
 		print OUT $currGene->toGFF()  if ( not defined( $gene_list{ $gene_id } ) );		
	}else{
 		print OUT $currGene->toGFF()  if $gene_list{ $gene_id } == 1;
	} 	
}

foreach my $curr_found ( keys %found ){
	print "NOT FOUND: $curr_found\n" if ( $found{$curr_found} ==  -1 );	
}
close(OUT);

print "Number of genes in the list: $num_genes_list\n";
print "Number of genes found: $num_genes_found\n";

	

