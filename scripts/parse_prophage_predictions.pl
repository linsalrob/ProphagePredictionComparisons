=pod

Parse the True Positive/True Negative output file, that I typically call *_tptn.tsv and summarize
the output as a table. The [python code](compare_predictions_to_phages.py) should just print this
out as a table, but it was repurporsed from something else and so it was easier to just write this
little parser.

Run with a directory of tptn files

=cut


use strict;
use Getopt::Std;
use Data::Dumper;
use FindBin qw($RealBin);
use lib "$RealBin/../ProphagePredictionsLib";
use Rob;
my $rob = new Rob;

=pod

Found:
Test set:
        Phage: 727 Not phage: 4898
Predictions:
        Phage: 844 Not phage: 4771
TP: 689  FP: 165  TN: 4733  FN: 38
Accuracy:    0.964      (this is the ratio of the correctly labeled phage genes to the whole pool of genes
Precision:   0.807      (This is the ratio of correctly labeled phage genes to all predictions)
Recall:      0.948      (This is the fraction of actual phage genes we got right)
Specificity: 0.966      (This is the fraction of non phage genes we got right)
f1 score: 0.872 (this is the harmonic mean of precision and recall, and is the best measure when, as in this case, there is a big difference between the number of phage and non-phage genes)

=cut

my %opts;
getopts('d:c:hv', \%opts);
unless ($opts{d} && $opts{c}) {
	die <<EOF;
	$0
	-d directory that contains the *tptn files (required)
	-c prophage caller (e.g. phispy, phigaro, virsorter, phage_finder (required)
	-h print header line
	-v verbose output


My standard naming convention is to put the sensitivity etc output in {genome}_{tool}_tptn.tsv

EOF
}


my @cols = ("TP", "TN", "FP", "FN", "Accuracy", "Precision", "Recall", "Specificity", "f1 score");

if ($opts{h}) {
	print join("\t", "Prophage Caller", "Genome", @cols), "\n";
}

opendir(DIR, $opts{d}) || die "can't open $opts{d}";
foreach my $f (grep {m/_tptn.tsv/} readdir(DIR)) {
	$f =~ /(\S+)_(\w+?)_tptn.tsv/;
	my $genome = $1;
	my $tool = $2;

	if ($f =~ /(\S+)_phage_finder_tptn.tsv/) {$genome = $1; $tool = "phage_finder"}
	if ($f =~ /(\S+)_phispy_trained_tptn.tsv/) {$genome = $1; $tool = "phispy_trained"}
	if ($f =~ /(\S+)_phispy_pvog_tptn.tsv/) {$genome = $1; $tool = "phispy_pvog"}

	unless (defined $genome) {print STDERR "Couldn't parse genome from $f\n"; die}
	unless ($tool eq $opts{c}) {print STDERR "You said you used $opts{c} but $opts{d} suggests it was $tool??\n"; die}

	my %res;
	map {$res{$_} = "NaN"} @cols;
	open(IN, "$opts{d}/$f") || die "can't open $opts{d}/$f";
	while (<IN>) {
		chomp;
		if (/^TP: (\d+)\s+FP: (\d+)\s+TN: (\d+)\s+FN: (\d+)/) {
			$res{"TP"} = $1;
			$res{"FP"} = $2;
			$res{"TN"} = $3;
			$res{"FN"} = $4;
		}
		if (/^(.*?):\s+([\d\.]+)/ && $res{$1}) {$res{$1}=$2}
	}
	close IN;
	print join("\t", "\t$tool", $genome, map{$res{$_}} @cols), "\n";
}

