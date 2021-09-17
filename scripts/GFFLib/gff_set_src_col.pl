#!/bin/env perl

# Set the 2nd column (src) of a GFF file to a value define by the user
use strict;

use FindBin;
use lib "$FindBin::Bin";

use GFFFile;


my $usage =
"gff_sort.pl. <src.  value> <GFF file> <output: modified GFF file>\n\n";

die $usage if scalar(@ARGV) != 3;

my $src_value		= $ARGV[0];
my $in_gff		= $ARGV[1];
my $out_gff		= $ARGV[2];


my $gffFile = GFFFile::new($in_gff);
$gffFile->read();

my $gffGenes = $gffFile->get_genes_hash();

# Ordering genes based on template name and start coord
my @gffGenesArray = values %{$gffGenes};


open OUT, ">$out_gff" or die "Unable to write on file $out_gff\n";
for my $currGene (@gffGenesArray) {
	my $content = $currGene->toGFF();
        my @lines = split "\n", $content;

	foreach my $curr_line ( @lines ){
		my @cols  = split "\t", $curr_line;
 		$cols[1] = $src_value;
		my $new_line = join "\t", @cols;
		print OUT $new_line . "\n";
	}
}
close(OUT);
