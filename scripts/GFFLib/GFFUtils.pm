package GFFUtils;
use strict 'vars';
use strict 'refs';

use GFFGene;
use GFFExon;
use GFFUTR;
use GFFCDS;
use GFFTranscript;

use Data::Dumper;

#use Set::IntervalTree;

sub are_end_exons_compatible {
	my ( $gffTrans1, $gffTrans2 ) = @_;

	# End exons are compatible if
	# - end coordinates of the first exon are the same
	# - start coordinates of the last exon are the same
	return 0
	  if $gffTrans1->get_first_exon()->get_end() !=
		  $gffTrans2->get_first_exon()->get_end();
	return 0
	  if $gffTrans1->get_end_exon()->get_start() !=
		  $gffTrans2->get_end_exon()->get_start();

	return 1;
}

sub are_exons_compatible {
	my ( $gffExon1, $gffExon2 ) = @_;

# Exons are compatible if:
# - They are the last exons and start coordinates are identical
# - They are the first exons and end coordinates are identical
# - If they are neither the first nor the last exons but the coordinates are identical

	return 1
	  if ( ( $gffExon1->get_start() == $gffExon2->get_start() )
		&& ( $gffExon1->get_end() == $gffExon2->get_end() ) );
	return 1
	  if ( $gffExon1->is_first()
		&& $gffExon2->is_first()
		&& ( $gffExon1->get_end() == $gffExon2->get_end() ) );
	return 1
	  if ( $gffExon1->is_last()
		&& $gffExon2->is_last()
		&& ( $gffExon1->get_start() == $gffExon2->get_start() ) );

	return 0;
}

sub ref_based_transcript_trimming {
	my ( $gffTransRef, $gffTransNew, $end ) =
	  @_;    #$end = [up,down], up = upstream end, down = downstream end

	if ( $end eq "up" ) {
		$gffTransNew->get_first_exon()
		  ->set_start( $gffTransRef->get_first_exon()->get_start() );
		$gffTransNew->get_first_UTR()
		  ->set_start( $gffTransRef->get_first_UTR()->get_start() );
		$gffTransNew->set_start( $gffTransRef->get_start() );
	}
	elsif ( $end eq "down" ) {
		my $new_exon_end = $gffTransRef->get_last_exon()->get_end();
		my $new_UTR_end  = $gffTransRef->get_last_UTR()->get_end();
		print "\tSetting to $new_exon_end - $new_UTR_end\n";
		$gffTransNew->get_last_exon()->set_end($new_exon_end);
		$gffTransNew->get_last_UTR()->set_end($new_UTR_end);
		$gffTransNew->set_end( $gffTransRef->get_end() );
	}
	else {
		die
"GFFUtils::ref_based_transcript_trimming third parameter should be either \"up\" or \"down\"\n";
	}

}

sub ref_based_exon_trimming {
	my ( $gffExonRef, $gffExonNew, $end ) =
	  @_;    #$end = [up,down], up = upstream end, down = downstream end

	return ref_based_exon_adjustment( $gffExonRef, $gffExonNew, $end );
}

sub ref_based_exon_adjustment {
	my ( $gffExonRef, $gffExonNew, $end ) = @_
	  ; #$end = [up,down], up = upstream end, down = downstream end, both = both ends will be adjusted (for single exons)

	my $gffTransNew = $gffExonNew->get_parent();
	my $gffTransRef = $gffExonRef->get_parent();

	if ( $end eq "up" || $end eq "both" ) {

		my $target_start = $gffExonRef->get_start();

		# Exon cannot shrink to a start value higher than the start attached CDS
		if ( $gffExonNew->has_CDS() ) {
			my $max_start = $gffExonNew->get_CDS()->get_start();
			$target_start = $max_start if ( $target_start > $max_start );
		}

		print ">> $target_start\n";
		$gffExonNew->set_start($target_start);
		$gffTransNew->set_start($target_start);

	}
	if ( $end eq "down" || $end eq "both" ) {

		my $target_end = $gffExonRef->get_end();

		# Exon cannot shrink to an end value lower than the end attached CDS
		if ( $gffExonNew->has_CDS() ) {
			my $min_end = $gffExonNew->get_CDS()->get_end();
			$target_end = $min_end if ( $target_end < $min_end );
		}

		print ">> $target_end\n";
		$gffExonNew->set_end($target_end);
		$gffTransNew->set_end($target_end);

	}

	if ( $end ne "up" && $end ne "down" && $end ne "both" ) {
		die
"GFFUtils::ref_based_exon_adjustment third parameter should be either \"up\" or \"down\"\n";
	}

}

#sub intersection{
#  my ($gffFILE1, $gffFILE2) = @_;
#
#  my $tree = Set::IntervalTree->new;
#
#  my $gffGenes1 = $gffFILE1->get_genes_hash();
#
#  for my $currGene (values %{$gffGenes1} ){
#	for my $currTranscript ( values %{$currGene->get_transcripts_hash()} ){
#		my $min = $currTranscript->get_exon()->get_start();
#		my $max = $currTranscript->get_exon()->get_end();
#
#		$tree->insert( new Number::Interval( Min=> $min, Max => $max ), $min, $max );
#
#	}
#  }
#
#  my $gffGenes2 = $gffFILE2->get_genes_hash();
#
#  for my $currGene (values %{$gffGenes} ){
#	for my $currTranscript ( values %{$currGene->get_transcripts_hash()} ){
#		my $min = $currTranscript->get_exon()->get_start();
#		my $max = $currTranscript->get_exon()->get_end();
#
#		my $ref_current_interval = new Number::Interval( Min=> $min, Max => $max );
#
#		# Retrieve intervals that overlap with current interval
#		my $arr_ref_intervals = $tree->fetch( $min, $max );
#
#		# Generate intersection interval
#		my $arr_ref_intersection_intervals = interval_intersection( $ref_current_interval, $arr_ref_intervals );
#
#	}
#  }
#
#}

sub _interval_intersection {
	my ( $ref_interval, $ref_arr_interval ) = @_;

	foreach my $curr_interval ( @{$ref_arr_interval} ) {
	}
}

sub sort_gene_arrays {
	my ( $gffGeneArrayRef, $numerical_chrom ) = @_;

	# Try to use the numerical part of chrom and then start coord of gene

	if ($numerical_chrom) {
		@{$gffGeneArrayRef} =
		  sort {

			my ($a_num_chrom) = ( $a->get_chrom() =~ /(\d+)/ );
			my ($b_num_chrom) = ( $b->get_chrom() =~ /(\d+)/ );

			return $a_num_chrom <=> $b_num_chrom
			  || $a->get_start() <=> $b->get_start()
		  } @{$gffGeneArrayRef};
	}
	else {
		@{$gffGeneArrayRef} =
		  sort {
			return $a->get_chrom() cmp $b->get_chrom()
			  || $a->get_start() <=> $b->get_start()
		  } @{$gffGeneArrayRef};

	}

}

sub is_gene_different_transcript_span {
	my ( $gene) = @_;

	if( $gene->get_transcript_span_start()  != $gene->get_start() ){
		return 1;
	}

	if( $gene->get_transcript_span_end()  != $gene->get_end() ){
		return 1;
	}

	return 0;
}

sub print_exon_comparison {
	my ( $gene_a, $gene_b ) = @_;

	my $gffTranscripts_a = $gene_a->get_transcripts_hash();
	my $gffTranscripts_b = $gene_b->get_transcripts_hash();

	for my $currTranscript_a ( values %{$gffTranscripts_a} ) {
		for my $currTranscript_b ( values %{$gffTranscripts_b} ) {
			for my $currExon_a ( @{ $currTranscript_a->get_exon_array() } ) {
				for my $currExon_b ( @{ $currTranscript_b->get_exon_array() } )
				{

					my $equivalence_symbol = "";
					if ( $currExon_a->identical($currExon_b) ) {
						$equivalence_symbol = "=";
					}
					elsif ( $currExon_a->overlaps($currExon_b) ) {
						$equivalence_symbol = "~";
					}
					print "\tEXON:\t"
					  . $currExon_a->get_id() . ":"
					  . $currExon_a->get_start() . "-"
					  . $currExon_a->get_end() . "\t"
					  . $equivalence_symbol . "\t"
					  . $currExon_b->get_id() . ":"
					  . $currExon_b->get_start() . "-"
					  . $currExon_b->get_end() . "\n"
					  if $equivalence_symbol ne "";
				}
			}
		}
	}
}

sub exon_comparison {
	my ( $gene_a, $gene_b ) = @_;

	my $content = "";
	
	my $gffTranscripts_a = $gene_a->get_transcripts_hash();
	my $gffTranscripts_b = $gene_b->get_transcripts_hash();

	for my $currTranscript_a ( values %{$gffTranscripts_a} ) {
		for my $currTranscript_b ( values %{$gffTranscripts_b} ) {
			for my $currExon_a ( @{ $currTranscript_a->get_exon_array() } ) {
				my $identical_or_overlaps = 0;

				for my $currExon_b ( @{ $currTranscript_b->get_exon_array() } )
				{
					$identical_or_overlaps = 1
					  if ( $currExon_a->identical($currExon_b)
						|| $currExon_a->overlaps($currExon_b) );
				}
				
				if ( $identical_or_overlaps == 0 ) {
					$content .= "\tEXON:\t"
					  . $currExon_a->get_id() . ":"
					  . $currExon_a->get_start() . "-"
					  . $currExon_a->get_end() . "\t" . "ONLY" . "\t" . "\n";
				}
			}
		}
	}

	for my $currTranscript_b ( values %{$gffTranscripts_b} ) {
		for my $currTranscript_a ( values %{$gffTranscripts_a} ) {
			for my $currExon_b ( @{ $currTranscript_b->get_exon_array() } ) {
				my $identical_or_overlaps = 0;
				for my $currExon_a ( @{ $currTranscript_a->get_exon_array() } )
				{
					$identical_or_overlaps = 1
					  if ( $currExon_b->identical($currExon_a)
						|| $currExon_b->overlaps($currExon_a) );
				}
				if ( $identical_or_overlaps == 0 ) {
					$content .= "\tEXON:\t\tONLY\t"
					  . $currExon_b->get_id() . ":"
					  . $currExon_b->get_start() . "-"
					  . $currExon_b->get_end() . "\n";
				}
			}
		}
	}

	for my $currTranscript_a ( values %{$gffTranscripts_a} ) {
		for my $currTranscript_b ( values %{$gffTranscripts_b} ) {
			for my $currExon_a ( @{ $currTranscript_a->get_exon_array() } ) {
				for my $currExon_b ( @{ $currTranscript_b->get_exon_array() } )
				{

					if ( not $currExon_a->identical($currExon_b)
						&& $currExon_a->overlaps($currExon_b) )
					{
						$content .= "\tEXON:\t"
						  . $currExon_a->get_id() . ":"
						  . $currExon_a->get_start() . "-"
						  . $currExon_a->get_end() . "\t" . "~" . "\t"
						  . $currExon_b->get_id() . ":"
						  . $currExon_b->get_start() . "-"
						  . $currExon_b->get_end() . "\n";
					}

				}
			}
		}
	}
	
	return $content;
}

sub overlapping_genes_ {
	my ( $gff_a, $gff_b ) = @_;

	return overlapping_genes_using_buffer( $gff_a, $gff_b, 0 );
}

sub overlapping_genes_using_buffer {

	# Reference GFFFile
	my ( $gff_a, $gff_b, $buffer ) = @_;

	my $gffAGenes = $gff_a->get_genes_hash();
	my $gffBGenes = $gff_b->get_genes_hash();

	my $gff_filename_a = $gff_a->get_filename();
	my $gff_filename_b = $gff_b->get_filename();

	# Ordering genes based on template name and start coord
	my @gffGenesArray = ( values %{$gffAGenes}, values %{$gffBGenes} );

	GFFUtils::sort_gene_arrays( \@gffGenesArray, 0 );

	my %overlap_from_A;
	my %overlap_from_B;

	for ( my $i1 = 0 ; $i1 < scalar(@gffGenesArray) ; $i1++ ) {

		#print $gffGenesArray[ $i1 ]->toGFF();
		#getc();

		for ( my $i2 = $i1 + 1 ; $i2 < scalar(@gffGenesArray) ; $i2++ ) {
			last
			  if (
				not $gffGenesArray[$i1]
				->overlaps( $gffGenesArray[$i2], "strandness", $buffer ) );
			next
			  if ( $gffGenesArray[$i1]->get_filename() eq
				$gffGenesArray[$i2]->get_filename() );

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

			#			print STDERR "GENE:\t"
			#			  . $gene_a->get_id() . ":"
			#			  . $gene_a->get_chrom() . ":"
			#			  . $gene_a->get_start() . "-"
			#			  . $gene_a->get_end() . ":"
			#			  . $gene_a->get_strand()
			#			  . ":exons="
			#			  . $gene_a->num_exons() . "\t"
			#			  . $gene_b->get_id() . ":"
			#			  . $gene_b->get_chrom() . ":"
			#			  . $gene_b->get_start() . "-"
			#			  . $gene_b->get_end() . ":"
			#			  . $gene_a->get_strand()
			#			  . ":exons="
			#			  . $gene_b->num_exons() . "\n";

			my $a_id = $gene_a->get_id();
			my $b_id = $gene_b->get_id();

			$overlap_from_A{$a_id} = $b_id;
			$overlap_from_B{$b_id} = $a_id;
		}

	}

	return ( \%overlap_from_A, \%overlap_from_B );
}

# List issues like this

# Current model                                       |||||-----||||||||----------||||||
# Final model                                         |||||-----||||||||||||||||||||||||

# OR

# Current model                                       |||||-----||||||||----------||||||
# Final model                                         |||||-----|||||||||||||||||||||

sub was_an_intron_retained {

	my ( $transA, $transB ) = @_;

	# Format:
	# <intron1 start>..<intron1 end>-<intron2 start>..<intron2 end>
	my $splice_str = $transB->get_splice_str();

	my @splices = split "-", $splice_str;

	foreach my $curr_splice (@splices) {
		my ( $splice_start, $splice_end ) =
		  ( $curr_splice =~ /(\d+)\.\.(\d+)/ );

		my $exon_A_overlaping_start = $transA->get_exon_by_coord($splice_start);
		my $exon_A_overlaping_end   = $transA->get_exon_by_coord($splice_end);

#print "$exon_A_overlaping_start $exon_A_overlaping_end\n" if ( $transcript_id eq 'AN7351_mRNA' );

		# If one of the attempts of retrieving the gene fails THEN
		# coord. still belongs to an intron. NO INTEGRAL intron retention
		next
		  if ( $exon_A_overlaping_start == 0 || $exon_A_overlaping_end == 0 );

		# If the same exon was retrieved using the $splice_start and $splice_end
		# THEN INTEGRAL RETENTION happened
		if ( $exon_A_overlaping_start->get_start() ==
			   $exon_A_overlaping_end->get_start()
			&& $exon_A_overlaping_start->get_end() ==
			$exon_A_overlaping_end->get_end() )
		{
			return 1;

		}
	}
	return 0;
}

# List issues like this

# Current model                                       |||||-----||||||||----------||||||
# Transcript evidence with partial intron retention   |||||-----|||||||||||
# Final model                                         |||||-----|||||||||||

sub was_an_intron_partially_retained {

	my ( $transA, $transB ) = @_;

# if none of the ends of current model is intronic when compared to previous version then
# this is not the case that we are looking for.
	if (   $transB->is_intronic( $transA->get_start() )
		|| $transB->is_intronic( $transA->get_end() ) )
	{
		return 1;
	}
	return 0;
}

return 1;

