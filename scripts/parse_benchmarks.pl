=pod

Parse out the benchmark data for the prophage callers and make a table

Rob Edwards, 2020

=cut

use strict;
use Getopt::Std;
use Data::Dumper;
use lib "$RealBin/../ProphagePredictionsLib";
use Rob;
my $rob = new Rob;

my %opts;
getopts('d:c:hv', \%opts);
unless ($opts{d} && $opts{c}) {
	die <<EOF;
	$0
	-d directory of benchmark results (required)
	-c prophage caller (required)
	-h print header
	-v verbose output
EOF
}

opendir(DIR, $opts{d}) || die "CAn't open $opts{d}";
foreach my $f (grep {$_ !~ /^\./} readdir(DIR)) {
	open(IN, "$opts{d}/$f") || die "can't open $opts{d}/$f";
	my $first  = 1;
	while (<IN>) {
		if ($first && $_ !~ /mean_load/) {
			print STDERR "$f does not appear to be a benchmark file. Skipped\n";
			close IN;
			next;
		}
		if ($first) {
			if ($opts{h}) {
				chomp;
				my @p = split /\s+/;
				print join("\t", "Prophage Caller", "Genome", @p), "\n";
				undef $opts{h};
			}
			$first = 0;
			next}
		chomp;
		my @p = split /\s+/;
		unless ($#p == 8) {
			print STDERR "$f does not appear to be a benchmark file. Too many columns in $_. Skipped\n";
			next;
		}
		print join("\t", $opts{c}, $f, @p), "\n";
	}
	close IN;
}





