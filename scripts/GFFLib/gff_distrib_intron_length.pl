#!/bin/env perl
use strict;

use FindBin;
use lib $FindBin::Bin;

use GFFFile;
use GFFUtils;
use Bio::SeqIO;
use Bio::Seq;

my $usage = "\nusage: $0  <GFF file> ";

my $debug = 0;

die $usage if ( scalar(@ARGV) != 1 );

my $gff_filename = $ARGV[0];

my $gffFile = GFFFile::new($gff_filename);

print "Reading GFF file...\n";
$gffFile->read( 1, 1 );

my $gffGenes = $gffFile->get_genes_hash();

# Ordering genes based on template name and start coord
my @gffGenesArray = values %{$gffGenes};
GFFUtils::sort_gene_arrays( \@gffGenesArray, 0 );

print "Number of genes in the GFF file: " . scalar(@gffGenesArray) . "\n";

my $cont_introns = 0;

my @introns_length;

for my $currGene (@gffGenesArray) {
	my $gene_id = $currGene->get_id();
	my $chrom   = $currGene->get_chrom();
	my $strand  = $currGene->get_strand();

	my $gffTranscripts = $currGene->get_transcripts_hash();
	for my $currTranscript ( values %{$gffTranscripts} ) {
		my $transcript_id = $currTranscript->get_id();

		my $intron_str = $currTranscript->get_splice_str();

		next if $intron_str eq "";

		my @introns = split "-", $intron_str;

		my $num_introns_gene = scalar(@introns);

		if ($debug) {
			print "Intron str: $intron_str Number of introns: "
			  . $num_introns_gene . "\n";
			getc();
		}

		for ( my $i = 0 ; $i < $num_introns_gene ; $i++ ) {
			my $currIntron = @introns[$i];

			my ( $start, $end ) = ( $currIntron =~ /(\d+)\.\.(\d+)/ );

			my $length = $end - $start + 1;

			push @introns_length, $length;
			if ($debug) {
				print
				  "$gene_id $transcript_id $cont_introns i= $i $start-$end\n";
				getc();
			}

			# Writing to files
			$cont_introns++;
		}

	}
}

# Calculate distribution

my @sorted_introns = sort { $a <=> $b } @introns_length;

my $ten_perc          = int( $cont_introns / 100 * 10 );
my $twenty_five_perc  = int( $cont_introns / 100 * 25 );
my $fifty_perc        = int( $cont_introns / 100 * 50 );
my $seventy_five_perc = int( $cont_introns / 100 * 75 );
my $ninety_perc       = int( $cont_introns / 100 * 90 );

my $num_sorted_introns = scalar(@sorted_introns);
if ($debug) {
	print "10% index:  $ten_perc  \n";
	print "25% index:  $twenty_five_perc  \n";
	print "50% index:  $fifty_perc  \n";
	print "75% index:  $seventy_five_perc \n";
	print "90% index:  $ninety_perc \n";

	print "Number of items: $num_sorted_introns\n";
}

print "Min:\t" . $sorted_introns[0] . "\n";
print "10%:\t" . $sorted_introns[$ten_perc] . "\n";
print "25%(Q1):\t" . $sorted_introns[$twenty_five_perc] . "\n";
print "50%:\t" . $sorted_introns[$fifty_perc] . "\n";
print "75%(Q3):\t" . $sorted_introns[$seventy_five_perc] . "\n";
print "90%:\t" . $sorted_introns[$ninety_perc] . "\n";
print "Max:\t" . $sorted_introns[ $num_sorted_introns - 1 ] . "\n";

my $q1 =  $sorted_introns[$twenty_five_perc];
my $q3 =  $sorted_introns[$seventy_five_perc];


# Lower bound: Q1 - 1.5 (Q1 - Q3)
my $outliers_lower_bound = $q1 - ( 1.5 * (  $q3 -  $q1 ) );

# Higher bound: Q3 + 1.5 (Q1 - Q3)
my $outliers_higher_bound = $q3 + ( 1.5 * (  $q3 -  $q1 ) );

$outliers_lower_bound = 0 if $outliers_lower_bound < 0; 

print "Lower bound:\t" . $outliers_lower_bound  . "\n";
print "Higher bound:\t" . $outliers_higher_bound  . "\n";
