use strict;
use Getopt::Std;
use Data::Dumper;
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
getopts('d:c:v', \%opts);
unless ($opts{d} && $opts{c}) {
	die <<EOF;
	$0
	-d directory that contains the *tptn files (required)
	-c prophage caller (e.g. phispy, phigaro, virsorter, phage_finder (required)
	-v verbose output


My standard naming convention is to put the sensitivity etc output in {genome}_{tool}_tptn.tsv

EOF
}


my @cols = ("Accuracy", "Precision", "Recall", "Specificity", "f1 score");

print join("\t", "Prophage Caller", "Genome", @cols), "\n";

opendir(DIR, $opts{d}) || die "can't open $opts{d}";
foreach my $f (grep {m/_tptn.tsv/} readdir(DIR)) {
	$f =~ /(\S+)_(\w+?)_tptn.tsv/;
	my $genome = $1;
	my $tool = $2;

	if ($f =~ /(\S+)_phage_finder_tptn.tsv/) {$genome = $1; $tool = "phage_finder"}

	unless (defined $genome) {print STDERR "Couldn't parse genome from $f\n"; die}
	unless ($tool eq $opts{c}) {print STDERR "You said you used $opts{c} but $opts{d} suggests it was $tool??\n"; die}

	my %res;
	map {$res{$_} = "NaN"} @cols;
	open(IN, "$opts{d}/$f") || die "can't open $opts{d}/$f";
	while (<IN>) {
		chomp;
		if (/^(.*?):\s+([\d\.]+)/ && $res{$1}) {$res{$1}=$2}
	}
	close IN;
	print join("\t", $tool, $genome, map{$res{$_}} @cols), "\n";
}

