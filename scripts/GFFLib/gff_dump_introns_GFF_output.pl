#!/bin/env perl
use strict;

use FindBin;
use lib "$FindBin::Bin";

use GFFFile;
use GFFUtils;
use Bio::SeqIO;
use Bio::Seq;

my $usage =
"\nusage: $0  <GFF file> <min. intron length> <out intron FASTA>";


die $usage if ( scalar(@ARGV) != 3 );

my $gff_filename               = $ARGV[0];
my $min_intron_length          = $ARGV[1];
my $out_intron                 = $ARGV[2];

my $gffFile = GFFFile::new($gff_filename);

print "Reading GFF file...\n";
$gffFile->read( 1, 1 );


my $gffGenes = $gffFile->get_genes_hash();

# Ordering genes based on template name and start coord
my @gffGenesArray = values %{$gffGenes};
GFFUtils::sort_gene_arrays( \@gffGenesArray, 0 );

print "Number of genes in the GFF file: " . scalar(@gffGenesArray) . "\n";

open OUT_INTRON, ">$out_intron" or die "Unable to open file $out_intron to write\n";

my $cont_introns   = 0;

for my $currGene (@gffGenesArray) {
	my $gene_id = $currGene->get_id();
	my $chrom   = $currGene->get_chrom();
	my $strand  = $currGene->get_strand();
	
	
	my $gffTranscripts = $currGene->get_transcripts_hash();
	for my $currTranscript ( values %{$gffTranscripts} ) {
		my $transcript_id  = $currTranscript->get_id();
		
		
		my $intron_str =  $currTranscript->get_splice_str();
		
		next if $intron_str eq "";
		
		my @introns = split "-", $intron_str;
		
		foreach my $currIntron ( @introns ){
			my ($start, $end) = ( $currIntron =~ /(\d+)\.\.(\d+)/ );
			
			my $length = $end - $start + 1;
			
			if( $length > $min_intron_length ){
				
								
				# Writing to files
				print OUT_INTRON "$chrom\tvoid\tintron\t$start\t$end\t.\t$strand\t.\tParent=$transcript_id\n";
				$cont_introns++;
			}
		}
		
	}
}

close(OUT_INTRON);
print "Number of intron sequences: $cont_introns\n";

