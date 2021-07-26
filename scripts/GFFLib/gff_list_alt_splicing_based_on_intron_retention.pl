#!/bin/env perl
use strict;

use FindBin;
use lib "$FindBin::Bin";

use GFFFile;
use GFFUtils;

# List issues like this

# Current model                                       |||||-----||||||||----------||||||
# Final model                                         |||||-----||||||||||||||||||||||||

# OR

# Current model                                       |||||-----||||||||----------||||||
# Final model                                         |||||-----|||||||||||||||||||||



# OR partial intron retentions


# Current model                                       |||||-----||||||||----------||||||
# Transcript evidence with partial intron retention   |||||-----|||||||||||
# Final model                                         |||||-----|||||||||||


my $usage =
"$0 <GFF file> <log file>\n";

die $usage if scalar(@ARGV) != 2;

my $gffFile = $ARGV[0];
my $log_file   = $ARGV[1];

print STDERR "Reading $gffFile ...\n";
my $gff_a = GFFFile::new($gffFile);
$gff_a->read();

my $gffAGenes = $gff_a->get_genes_hash();

open LOGFILE, ">$log_file" or die "Unable to open logfile $log_file\n";
print STDERR "Checking alt. splicing ...\n";
for my $gene_id ( keys %{$gffAGenes} ) {

	my $gffTranscriptsA = $gffAGenes->{$gene_id}->get_transcripts_hash();
		
	my $mainTranscript;
	map {$mainTranscript =  $gffTranscriptsA->{ $_ } if $_ =~ /.1$/ } keys %{$gffTranscriptsA};
	
	for my $transcript_id ( keys %{$gffTranscriptsA} ) {
		next if $transcript_id eq $mainTranscript->get_id();
		
		my $altSplicing = $gffTranscriptsA->{$transcript_id};
		
		if( GFFUtils::was_an_intron_partially_retained( $altSplicing, $mainTranscript  ) ){
				print LOGFILE $altSplicing->get_id() . "\n";
				print OUT $altSplicing->get_id() . "\n";
				next;			
		}
		
		if( GFFUtils::was_an_intron_retained( $altSplicing, $mainTranscript  ) ){
				print LOGFILE $altSplicing->get_id() . "\n";
				print OUT $altSplicing->get_id() . "\n";
				next;			
		}		
	}
}

close(LOGFILE);
