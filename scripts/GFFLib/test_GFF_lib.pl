#!/usr/bin/perl
use strict;

use FindBin;
use lib $FindBin::Bin;

use GFFFile;
use GFFUtils;

my $usage = "test_GFF_lib.pl <GFF A>\n\n";

die $usage if scalar(@ARGV) != 1;

my $gff_filename_a = $ARGV[0];

my $gff_a = GFFFile::new($gff_filename_a);
$gff_a->read();

#################
# Print chromosomes in the GFF

foreach my $curr_chrom ( $gff_a->get_chrom_names() ){
	print $curr_chrom . "\n";	
}