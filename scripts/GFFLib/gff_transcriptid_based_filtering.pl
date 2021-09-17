#!/bin/env perl
use strict;

use FindBin;
use lib "$FindBin::Bin";

use GFFFile;

my $usage = "$0 <id list file> <GFF file> <new GFF> [remove]\n\n";

die $usage if ( scalar(@ARGV) != 3 && scalar(@ARGV) != 4 );

my $idFile = $ARGV[0];
my $inFile = $ARGV[1];
my $outFile = $ARGV[2];
my $shouldIremove = $ARGV[3];


# Indicate if the transcript list should be either kept or removed from the GFF file
my $remove = 0;
$remove = 1 if ( defined( $shouldIremove ) );

# Transcripts list
my %transcript_list;

my $num_transcripts_list = 0;
my %found;

# Read file listing transcripts
open IN, $idFile or die "Unable to open file $idFile\n";
while (<IN>) {
	my $line = $_;
	chomp $line;
	$transcript_list{$line} = 1;
	$num_transcripts_list++;
	$found{$line} = -1;
}

# Reading GFF file
my $gffFile = GFFFile::new( $inFile );
$gffFile->read();

my $gffGenes = $gffFile->get_genes_hash();

my $num_transcripts_found = 0;


open OUT, ">$outFile" or die "Unable to write on file $outFile\n";
for my $currGene ( values %{$gffGenes} ) {
	for my $currTranscript ( values %{ $currGene->get_transcripts_hash() } ) {

		my $transcript_id = $currTranscript->get_id();

		if ( $transcript_list{$transcript_id} == 1 ) {
			$found{$transcript_id} = 1;
			$num_transcripts_found++;

			$currGene->delete_transcript($transcript_id);
						
		}
	}
	print OUT $currGene->toGFF();	
}

foreach my $curr_found ( keys %found ) {
	print STDERR "NOT FOUND: $curr_found\n" if ( $found{$curr_found} == -1 );
}

print STDERR "Number of transcripts in the list: $num_transcripts_list\n";
print STDERR "Number of transcripts found: $num_transcripts_found\n";

close(OUT);
