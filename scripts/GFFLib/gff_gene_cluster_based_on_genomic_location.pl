#!/bin/env perl
use strict;


use FindBin;
use lib $FindBin::Bin;


use GFFFile;
use GFFUtils;



my $debug = 0;

# Define cluster of genes based on their location in the genome (gene locus).
# Genes that overlap are considered members of the same cluster
# The gene with the longest CDS is chosen as the representative of the cluster
# All representatives are reported in the GFF OUT file


my $usage = "$0 <GFF IN> <GFF OUT with one representative for each cluster > \n\n";

die $usage if scalar(@ARGV) != 2;

my $gff_filename_in = $ARGV[0];
my $out_gff         = $ARGV[1];

my $gff_handler = GFFFile::new( $gff_filename_in );
$gff_handler->read();


my $gffGenes = $gff_handler->get_genes_hash();

# Ordering genes based on template name and start coord
my @gffGenesArray = values %{$gffGenes};
GFFUtils::sort_gene_arrays( \@gffGenesArray, 0 );


# As genes get added to the cluster the program adjust those variables
# They represent left and right end of the range spanning all genes in
# the clusters

my $cluster_start = -1;
my $cluster_end   = -1;

my $last_gene_chrom;

my @cluster;

open OUT, ">$out_gff" or die "Unable to open file $out_gff to write\n";

for my $currGene (@gffGenesArray) {
	my $gene_id        = $currGene->get_id();
	my $gffTranscripts = $currGene->get_transcripts_hash();
	my $start =  $currGene->get_start();
	my $end   =  $currGene->get_end();
	my $chrom =  $currGene->get_chrom();
	
	
	print STDERR "$gene_id...\n" if $debug;
	
	if( $start > $end ){
		my $temp = $start;
		$start = $end;
		$end = $temp;
	}
	
	
	# Start new clusters if no intersection
	# or if chromosome have changed
	if ( $start > $cluster_end || $chrom ne $last_gene_chrom ){
				
		# If previous cluster not empty, choose the representative
		# gene and dump it
		if( scalar ( @cluster ) != 0 ){
			my $chosen_gene = choose_representative( \@cluster );
			print OUT $chosen_gene->toGFF();
			undef @cluster;
			$cluster_start = -1;
			$cluster_end   = -1;
			
		}		   	   		   
		$cluster_start = $start;
		
	}

	# Extend end of the cluster if needed
	$cluster_end   = $end if $end > $cluster_end;
	push @cluster, $currGene;
	
	
	$last_gene_chrom = $chrom;

}
my $chosen_gene = choose_representative( \@cluster );
print OUT $chosen_gene->toGFF();


close(OUT);



exit(0);


sub choose_representative{
	my ($clusterRef) = @_;
	
	my $chosen;
	my $longest_cds = 0;
	
	foreach my $currGene ( @{$clusterRef} ){
		my $id = $currGene->get_id();

		my $cds_length = 0;

		my $gffTranscripts = $currGene->get_transcripts_hash();
		
		
		for my $currTranscript ( values %{$gffTranscripts} ) {		
			my $temp_cds_length = $currTranscript->get_CDS_length();
			$cds_length = $temp_cds_length if $temp_cds_length > $cds_length;
		}
		
		print STDERR "ID: $id  CDS_LENGTH: $cds_length\n";
		
		if( $cds_length > $longest_cds ){
			$chosen = $currGene;
			$longest_cds = $cds_length;
		}
			
		
	}
	my $chosen_id = $chosen->get_id();
	print STDERR "CHOSEN ID: $chosen_id\n";
	print STDERR "======================\n";
	
	return $chosen;
}




