#!/usr/bin/perl
use strict;

use FindBin;
use lib $FindBin::Bin;

use GFFFile;
use GFFUtils;

my $usage = "gff_compare.pl <GFF A> <GFF B>\n\n";

die $usage if scalar(@ARGV) != 2;

my $gff_filename_a = $ARGV[0];
my $gff_filename_b = $ARGV[1];

my $gff_a = GFFFile::new($gff_filename_a);
$gff_a->read();

my $gff_b = GFFFile::new($gff_filename_b);
$gff_b->read();

my $gffAGenes = $gff_a->get_genes_hash();
my $gffBGenes = $gff_b->get_genes_hash();

# Ordering genes based on template name and start coord
my @gffGenesArray = ( values %{$gffAGenes}, values %{$gffBGenes} );

GFFUtils::sort_gene_arrays( \@gffGenesArray, 0 );

# This hash prevent a second comparison of genes that have identical
# counterpart in the other file
my %identicals;

# This check will keep track of genes that at least overlaps with other
# genes. Later it will be used to report genes that are only found
# in each of the files
my %similars;

for ( my $i1 = 0 ; $i1 < scalar(@gffGenesArray) ; $i1++ ) {

	#print $gffGenesArray[ $i1 ]->toGFF();
	#getc();

	next if defined $identicals{ $gffGenesArray[$i1]->get_id() };

	for ( my $i2 = $i1 + 1 ; $i2 < scalar(@gffGenesArray) ; $i2++ ) {
		last if ( not $gffGenesArray[$i1]->overlaps( $gffGenesArray[$i2] ) );
		next
		  if ( $gffGenesArray[$i1]->get_filename() eq
			$gffGenesArray[$i2]->get_filename() );
		next if defined $identicals{ $gffGenesArray[$i2]->get_id() };

		my $gene_a;
		my $gene_b;

		if (   $gffGenesArray[$i1]->get_filename() eq $gff_filename_a
			&& $gffGenesArray[$i2]->get_filename() eq $gff_filename_b )
		{
			$gene_a = $gffGenesArray[$i1];
			$gene_b = $gffGenesArray[$i2];

		}
		elsif ($gffGenesArray[$i1]->get_filename() eq $gff_filename_b
			&& $gffGenesArray[$i2]->get_filename() eq $gff_filename_a )
		{
			$gene_b = $gffGenesArray[$i1];
			$gene_a = $gffGenesArray[$i2];

		}
		else {
			die "Unknown filenames:\n\t"
			  . $gffGenesArray[$i1]->get_filename() . "\n\t"
			  . $gffGenesArray[$i2]->get_filename() . "\n";
		}

		my $exon_comparison = GFFUtils::exon_comparison( $gene_a, $gene_b );

		if ( $exon_comparison ne '' ) {
			print "GENE:\t"
			  . $gene_a->get_id() . ":"
			  . $gene_a->get_chrom() . ":"
			  . $gene_a->get_start() . "-"
			  . $gene_a->get_end() . ":"
			  . $gene_a->get_strand()
			  . ":exons="
			  . $gene_a->num_exons() . "\t"
			  . $gene_b->get_id() . ":"
			  . $gene_b->get_chrom() . ":"
			  . $gene_b->get_start() . "-"
			  . $gene_b->get_end() . ":"
			  . $gene_a->get_strand()
			  . ":exons="
			  . $gene_b->num_exons() . "\n";
			print $exon_comparison;
			#print "SIM\t" . $gene_a->get_id() . "\n";
			#print "SIM\t" . $gene_b->get_id() . "\n";
			$similars{ $gene_a->get_id() } = 1;
			$similars{ $gene_b->get_id() } = 1;
		}
		else {
			#print "ID\t" . $gene_a->get_id() . "\n";
			#print "ID\t" . $gene_b->get_id() . "\n";
			$identicals{ $gene_a->get_id() } = 1;
			$identicals{ $gene_b->get_id() } = 1;
			last;
		}

	}

}

# Report genes only found in one of the GFFs

for ( my $i1 = 0 ; $i1 < scalar(@gffGenesArray) ; $i1++ ) {

	my $gene    = $gffGenesArray[$i1];
	my $gene_id = $gene->get_id();

	if (   ( not defined( $identicals{$gene_id} ) )
		and ( not defined( $similars{$gene_id}   ) )
	   )
	{
		#print "ONLY\t" . $gene_id . "\n";
		
		my $content =
		    $gene->get_id() . ":"
		  . $gene->get_chrom() . ":"
		  . $gene->get_start() . "-"
		  . $gene->get_end() . ":"
		  . $gene->get_strand()
		  . ":exons="
		  . $gene->num_exons();

		if ( $gene->get_filename() eq $gff_filename_a ) {
			print "GENE:\t" . $content . "\t" . "ONLY\n";
		}
		elsif ( $gene->get_filename() eq $gff_filename_b ) {
			print "GENE:\t" . "ONLY" . "\t" . $content . "\n";
		}
		else {
			die "Unknown filename:\n\t" . $gene->get_filename() . "\n";
		}
	}

}
