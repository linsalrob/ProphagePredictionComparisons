#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
my $h='';
my $in_dir='';
GetOptions ('help' => \$h, 'h' => \$h, 'i=s'=>\$in_dir);
if ($h==1 || $in_dir eq ""){ # If asked for help or did not set up any argument
	print "# Script to reverse-engineer the identifiers from VirSorter1, because VirSorter1 tends to change the contig names...
# Arguments :
# -i: input virsorter results folder
";
	die "\n";
}

my $translation_file=$in_dir."/fasta/input_sequences_id_translation.tsv";
my $locs_file=$in_dir."/locs.tsv.vsidentifiers";
my $out_file=$in_dir."/locs.tsv";

my %hash;
open my $tsv,"<",$translation_file;
while(<$tsv>){
      chomp($_);
      my @tab=split("\t",$_);
      $tab[1]=~s/^VIRSorter_//;
      $hash{$tab[1]}=$tab[0];
}
close $tsv;

open my $s1,">",$out_file;
open my $tsv,"<",$locs_file;
while(<$tsv>){
      chomp($_);
      my @tab=split("\t",$_);
      if (!defined($hash{$tab[0]})){print "Pblm - $tab[0] should have a translation from $translation_file\n"; die("\n");}
      $tab[0]=$hash{$tab[0]};
      print $s1 join("\t",@tab)."\n";
}
close $tsv;
close $s1;
