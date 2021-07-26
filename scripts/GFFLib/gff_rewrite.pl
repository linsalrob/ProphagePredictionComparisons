#!/usr/bin/env perl

# TODO:
# * Keep the comments in the beginning of the file (DONE)
# * Disregard FASTA sequences in the end of the file (DONE)
# * Implement getOpts (DONE)
# * Check if this functionality is working: --trim_genes_based_on_transcripts
# * Propagate Name attribute to gene record (from CDS record) (POSTPONED)
# * Leave unchanged the second column. (POSTPONED)
# * If keep_id, and missing mRNA or gene: propagate the ID to all other records. Double check if ID was not used (POSTPONED)

use strict;
use GFFFile;
use GFFUtils;
use Carp;
use Pod::Usage;
use Getopt::Long;

$Carp::Verbose = 1;

=head1 NAME

gff_rewrite.pl

=head1 SYNOPSIS

gff_rewrite.pl 
	--input <GFF IN> 
	--output <GFF OUT> 
	[  --change_id_single_exon_id | 
		--change_id_multi_exon_id ] 
	[ --discard_non_essential_attributes ]
	[ --add_missing_features ]		

=head1 OPTIONS

B<--input> - GFF input

B<--output> - GFF output

B<--change_id_single_exon_id> - B<(Optional)> Indicates that IDs should be normalized according to the following rules: 
	transcript id = <gene_id>-T;
	CDS id = <gene_id>-P;
	Exon id = <gene_id>-E;
	UTR id = <gene_id>-UTR
	
B<--change_id_multi_exon_id> - B<(Optional)> Same as change_id_single_exon_id, but in this case the exon IDs will 
be numbered (e.g. <gene_id>-E1, <gene_id>-E2)

B<--discard_non_essential_attributes> - B<(Optional)> Maintain only the fields Parent, ID and Name on column 9 

B<--add_missing_features> -  B<(Optional)> Add missing exons based on corresponding CDS, missing genes from 
transcript, missing transcripts from exons and missing exons from transcripts (only non-coding features)

B<--help> - prints the usage information. B<(Optional)>

=head1 DESCRIPTION


=head1 CONTACT
 Gustavo C. Cerqueira (2018)
 cerca11@gmail.com
 gcerqueira@pgdx.com
=cut

my ($in_gff, $out_gff, $change_id_single_exon_id, $change_id_multi_exon_id );
my ($discard_non_essential_attributes, $add_missing_features);
my $help;

GetOptions(	'input=s'							=> \$in_gff,
			'output=s'							=> \$out_gff,
			'change_id_single_exon_id'			=> \$change_id_single_exon_id,
			'change_id_multi_exon_id'			=> \$change_id_multi_exon_id,
			'discard_non_essential_attributes' 	=> \$discard_non_essential_attributes,
			'add_missing_features' 				=> \$add_missing_features,
			'help=i'								=> \$help );
		
		

if( not defined($in_gff) ){
	pod2usage(
		-message => "Error: Parameter --input is required!\n\n",
		-verbose => 1,
		-exitval => 1,
		-output  => \*STDERR
	);   
} 

if( not defined($out_gff) ){
	pod2usage(
		-message => "Error: Parameter --output is required!\n\n",
		-verbose => 1,
		-exitval => 1,
		-output  => \*STDERR
	);   
} 


if( defined($change_id_single_exon_id) && defined($change_id_multi_exon_id) ){
	pod2usage(
		-message => "Error: Either use --change_id_single_exon_id or --change_id_multi_exon_id !\n\n",
		-verbose => 1,
		-exitval => 1,
		-output  => \*STDERR
	);   
} 


if( defined($help) ){
   pod2usage(-verbose => 1 ,-exitval => 2);
} 


my $force_id_change = 0;
if( defined($change_id_single_exon_id) || defined($change_id_multi_exon_id) ){
	$force_id_change = 1;
}

my $gffFile = GFFFile::new($in_gff);
$gffFile->read($discard_non_essential_attributes,$add_missing_features);
my $gffGenes = $gffFile->get_genes_hash();

# Ordering genes based on template name and start coord
my @gffGenesArray = values %{$gffGenes};
GFFUtils::sort_gene_arrays( \@gffGenesArray, 0 );

open OUT, ">$out_gff" or croak "ERROR: Unable to open file $out_gff to write\n";

print OUT $gffFile->get_header_comments();

for my $currGene (@gffGenesArray) {
		my $gene_id        = $currGene->get_id();
		$currGene->set_name( $gene_id )  if ($force_id_change);
		my $gffTranscripts = $currGene->get_transcripts_hash();
		for my $currTranscript ( values %{$gffTranscripts} ) {

			$currTranscript->force_multi_cds();
			$currTranscript->correct_cds_phase();


			$currTranscript->set_id( $gene_id . "-T" ) if ($force_id_change);
			
			# Set feature type to mRNA if it was originally as "transcript"
			$currTranscript->set_transcript_type( "mRNA")  
			  if $currTranscript->get_transcript_type() eq "transcript";

			my $count = 1;
			for my $currExon ( @{ $currTranscript->get_exon_array() } ) {
				$currExon->set_id( $gene_id . "-E" . $count ) 
					if ( defined($change_id_multi_exon_id) );
					
				$currExon->set_id( $gene_id . "-E"  ) 
					if ( defined($change_id_single_exon_id) );
				$count++;
			}

			$count = 1;
			for my $currCDS ( @{ $currTranscript->get_CDS_array() } ) {
				$currCDS->set_id( $gene_id . "-P" ) if ($force_id_change);
				$count++;
			}

			$count = 1;
			for my $currUTR ( @{ $currTranscript->get_UTR_array() } ) {
				$currUTR->set_id( $gene_id . "-UTR" . $count ) if ($force_id_change);
				$count++;
			}
			
			$currTranscript->force_multi_cds();
			$currTranscript->correct_cds_phase();
		}
	
	print OUT $currGene->toGFF();	
}
close(OUT);

