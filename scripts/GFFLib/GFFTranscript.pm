package GFFTranscript;
use strict 'vars';
use strict 'refs';

use Carp;

use GFFExon;
use GFFUTR;
use GFFCDS;
#use Interval;

sub new {
	my ( $transcript_id, $transcript_type, $start_coord, $end_coord, $parent ) =
	  @_;

	my $self = {
		id              => $transcript_id,
		transcript_type => $transcript_type,
		start           => $start_coord,
		end             => $end_coord,
		parent          => $parent
	};

	bless $self, GFFTranscript;

	return $self;
}

sub get_parent {
	my $self = shift;
	return $self->{parent};
}

sub get_transcript_type {
	my $self = shift;
	return $self->{transcript_type};
}

sub set_transcript_type {
	my $self = shift;
	my ($new_type) = @_;
	$self->{transcript_type} = $new_type;
}

# Generate a string summnaryzing all splice junctions
# within a certain pair of coordinates
# Format:
# <intron1 start>..<intron1 end>-<intron2 start>..<intron2 end>
sub get_splice_str_within {
	my $self = shift;

	my ( $rng_start, $rng_end ) = @_;

	my $exon_arr  = $self->{exons};
	my $num_exons = scalar( @{$exon_arr} );

	return "" if $num_exons == 1;

	my $str = "";

	for ( my $exon_ind = 0 ; $exon_ind < ( $num_exons - 1 ) ; $exon_ind++ ) {

		my $intron_start = $exon_arr->[$exon_ind]->get_end() + 1;
		my $intron_end   = $exon_arr->[ $exon_ind + 1 ]->get_start() - 1;

		next if $intron_start < $rng_start;
		last if $intron_end > $rng_end;

		$str .= $intron_start . ".." . $intron_end . "-";
	}

	$str =~ s/-$//;

	return $str;

}

sub is_intronic {
	my $self = shift;

	my ($coord) = @_;

	my $exon_arr  = $self->{exons};
	my $num_exons = scalar( @{$exon_arr} );

	return 0 if $num_exons == 1;

	for ( my $exon_ind = 0 ; $exon_ind < ( $num_exons - 1 ) ; $exon_ind++ ) {

		my $intron_start = $exon_arr->[$exon_ind]->get_end() + 1;
		my $intron_end   = $exon_arr->[ $exon_ind + 1 ]->get_start() - 1;

		return 1 if ( $coord >= $intron_start && $coord <= $intron_end );
	}

	return 0;
}

# Generate a string summnaryzing all splice junctions
# Format:
# <intron1 start>..<intron1 end>-<intron2 start>..<intron2 end>
sub get_splice_str {
	my $self = shift;

	my $exon_arr  = $self->{exons};
	my $num_exons = scalar( @{$exon_arr} );

	return "" if $num_exons == 1;

	my $str = "";

	for ( my $exon_ind = 0 ; $exon_ind < ( $num_exons - 1 ) ; $exon_ind++ ) {

		my $intron_start = $exon_arr->[$exon_ind]->get_end() + 1;
		my $intron_end   = $exon_arr->[ $exon_ind + 1 ]->get_start() - 1;

		$str .= $intron_start . ".." . $intron_end . "-";
	}

	$str =~ s/-$//;

	return $str;

}

sub attach_UTR {
	my $self = shift;

	# Attach cds utr to exon
	foreach my $ref_utr ( @{ $self->{UTR} } ) {
		foreach my $ref_exon ( @{ $self->{exons} } ) {
			if (   $ref_exon->get_start() == $ref_utr->get_start()
				|| $ref_exon->get_end() == $ref_utr->get_end() )
			{
				#print "Attaching UTR " . $ref_utr->get_id() . " to exon " . $ref_exon->get_id() . "\n";  
				$ref_exon->attach_UTR($ref_utr);
				last;
			}
		}
	}
}

sub attach_CDS {
	my $self = shift;

	foreach my $ref_cds ( @{ $self->{CDS} } ) {
		foreach my $ref_exon ( @{ $self->{exons} } ) {
			if (   $ref_exon->get_start() <= $ref_cds->get_start()
				&& $ref_exon->get_end() >= $ref_cds->get_end() )
			{

				if ( $ref_exon->has_CDS() ) {
					print STDERR "Error at trying to attach CDS \'"
					  . $ref_cds->get_id()
					  . "\' to exon \'"
					  . $ref_exon->get_id()
					  . "\'. This exon is already attached to CDS: \'"
					  . $ref_exon->get_CDS()->get_id() . "\'\n";
					print STDERR "Possibly a duplicated gene!\n";
					print STDERR "Continuing processing of file!\n";
					next;
				}

				$ref_exon->attach_CDS($ref_cds);
				last;
			}
		}
	}

}

sub sort_exon_arr {
	my $self = shift;
	
	croak "ERROR: Transcript " . $self->get_id() . " does not have exons!!!" 
		if scalar( $self->{exons} ) == 0; 

	# Sort based on absolute coordinates
	@{ $self->{exons} } =
	  sort { return $a->get_start() <=> $b->get_start() } @{ $self->{exons} };

	# Sort  based on coding strand coordinates
	if ( $self->get_strand() eq "-" ) {
		@{ $self->{exons_strand} } = reverse @{ $self->{exons} };
	}
	else {
		$self->{exons_strand} = $self->{exons};
	}

}

sub sort_UTR_arr {
	my $self = shift;

	# Sort based on absolute coordinates
	@{ $self->{UTR} } =
	  sort { return $a->get_start() <=> $b->get_start() } @{ $self->{UTR} };

	# Sort  based on coding strand coordinates
	if ( $self->get_strand() eq "-" ) {
		@{ $self->{UTR_strand} } = reverse @{ $self->{UTR} };
	}
	else {
		$self->{UTR_strand} = $self->{UTR};
	}

}

sub sort_CDS_arr {
	my $self = shift;

	# Sort based on absolute coordinates
	@{ $self->{CDS} } =
	  sort { return $a->get_start() <=> $b->get_start() } @{ $self->{CDS} };

	# Sort  based on coding strand coordinates
	if ( $self->get_strand() eq "-" ) {
		@{ $self->{CDS_strand} } = reverse @{ $self->{CDS} };
	}
	else {
		$self->{CDS_strand} = $self->{CDS};
	}

}

sub post_processing {
	my $self = shift;
	
	my ($add_missing) = @_;
		
	$self->attach_UTR();
	$self->attach_CDS();

	 
	if( $add_missing == 1 ){
		$self->add_missing_exons_from_CDS()	
	}
		
	# If non-coding and it doesn't have an exon, then add one
	if( $self->num_exons() == 0  && $self->num_CDSs() == 0 && $add_missing == 1  ){
		$self->add_missing_exons_from_transcript();
	}
	
	#$self->add_UTR_element_based_on_CDS_exon_difference();


	# Sort exon, CDS and UTR based on absolute coordinates and strand based coordinates
	$self->sort_exon_arr();
	$self->sort_CDS_arr();
	$self->sort_UTR_arr();

	
}


sub add_missing_exons_from_CDS {
	my $self = shift;

	foreach my $ref_CDS ( @{ $self->{'CDS'} } ) {
		if( not $ref_CDS->is_attached() ){
			

			my $gene_id     = $self->get_parent()->get_id();
			my $exon_id     = $gene_id . "-E";
					

			my $ref_exon = $self->add_exon( $exon_id, $ref_CDS->get_start(), $ref_CDS->get_end() );
			
			$self->get_parent()->get_file()->add_exon_to_index( $gene_id, $self, $self->get_parent() ); 
			
			$ref_exon->attach_CDS($ref_CDS);
		}
		
	}
}


sub add_missing_exons_from_transcript {
	my $self = shift;

			

			my $gene_id     = $self->get_parent()->get_id();
			my $exon_id     = $gene_id . "-E";
					

			my $ref_exon = $self->add_exon( $exon_id, $self->get_start(), $self->get_end() );
			
			$self->get_parent()->get_file()->add_exon_to_index( $gene_id, $self, $self->get_parent() ); 
			
}


#sub add_UTR_element_based_on_CDS_exon_difference {
#	my $self = shift;
#	
#	# Return if its non-coding gene
#	return if $self->num_CDSs() == 0;
#		
#	for my $currExon ( @{$self->get_exon_array()} ){
#							
#			my $exon_id = $currExon->get_id();						
#			# If current exon does not have CDS neither a UTR attached to it
#			# this indicates a non-coding exon without an UTR attached
#			# ADD the UTR!
#			if( not $currExon->has_CDS() && not $currExon->has_UTR()  ){
#				
#				my $utr_type;
#				
#				my $cds_length_on_left = 
#				  $self->get_CDS_length_within( $self->get_start(), $currExon->get_start() - 1 );
#
#				my $cds_length_on_right = 
#				  $self->get_CDS_length_within( $currExon->get_end() + 1 $self->get_end() );
#				  
#				croak "ERROR: Trying to create an UTR element in the non-coding exon $exon_id\n" . 
#				      "This exon is non-coding but is flanked both upstream and downstream by coding exons!\n\n" 
#				if( $cds_length_on_left > 0 || $cds_length_on_rigth > 0 );
#				
#				
#				
#				if( $cds_length_on_right > 0 && $self->get_strand() eq '+' ||
#				    $cds_length_on_left > 0 && $self->get_strand() eq '-'
#				    ){
#					$utr_type = 'five_prime_utr'
#				 }else{
#				 	$utr_type = 'three_prime_utr'
#				 }
#				 
#				
#				# Check if all exons either upstream or downstream
#				# do NOT have a CDS. Otherwise issue an error.
#				if (  == 0 ||
#				      == 0 )
#				
#				# Create an UTR representing the whole exon
#				
#			
#			}
#			#If exon 5' end different than CDS boundaries
#			if{				
#				# Check if all upstream exons do not have CDS
#				# Otherwise issue an error
#				
#				
#				# Create a UTR from the start of the exon
#				# to the the start of the CDS - 1
#				
#				
#				
#			}
#			#If exon 3' end different than CDS boundaries
#			if{				
#				# Check if the differences is due to STOP codon
#				# disregard if that's the case
#
#				# Check if all downstream exons do not have CDS
#				# Otherwise issue an error
#				
#				
#				# Create a UTR from the start of CDS + 1 to
#				# the end of the exon
#				
#								
#				
#			}
#				 
#				
#				 
#			
#			my $currCDS = $currExon->get_CDS();
#			
#			if( $currCDS->{start} !=  $currExon->{start} || $currCDS->{end} != $currExon->{end} ){
#				
#				# Checking if the difference is due to STOP codon
#				next if( abs( $currCDS->{end} -  $currExon->{end} ) == 3 && $currExon->is_last() && $currGene->get_strand() eq "+" );
#				next if( abs( $currCDS->{start} -  $currExon->{start} ) == 3 && $currExon->is_first() && $currGene->get_strand() eq "-" );
#				
#				print "\tTranscript: " . $currTranscript->get_id() . 
#				      " Strand: " . $currGene->get_strand() . 
#				      " " . $currTranscript->get_exon_span_start() . 
#				      " " . $currTranscript->get_exon_span_end() . "\n";
#				      
#				print "\t\tCDS  Start:" . $currCDS->{start} . "\t\tEnd:" . $currCDS->{end} . "\n";
#				print "\t\tExon  Start:" . $currExon->{start} . "\t\tEnd:" . $currExon->{end} . "\n";
#				print "\n";
#				print $currGene->toGFF();
#				print "\n";
#				getc();
#			
#			}elsif( $currExon
#				
#				
#				
#				
#				elsif( $currTranscript->get_exon_span_start() != $currGene->get_start() ||
#			        $currTranscript->get_exon_span_end()   != $currGene->get_end()      ){
#			    print "Different>>>\n";
#			    print $currGene->toGFF();
#			    getc();
#			}
#				
#		}		 	
#	
#	
#	my $five_prime_exon = 
#	if( $self->
#
#	
#		if ( $feature eq $UTR5 ) {
#			my $curr_gene_name =
#			  $self->get_gene_name_by_transcript_name($parent);
#
#			# If it doesn't have a name create one
#			if ( $id eq "" ) {
#				$id =
#				  $parent . "-UTR"
#				  . ( $self->{genes}->{$curr_gene_name}->get_transcript($parent)
#					  ->num_UTRs() + 1 );
#			}
#			print STDERR
#			  "UTR5: $id\tTranscript: $parent\tGene: $curr_gene_name\n"
#			  if $debug;
#
#			$temp_feature =
#			  $self->{genes}->{$curr_gene_name}->get_transcript($parent)
#			  ->add_5UTR( $id, $start, $end );
#		}
#
#		if ( $feature eq $UTR3 ) {
#			my $curr_gene_name =
#			  $self->get_gene_name_by_transcript_name($parent);
#
#			# If it doesn't have a name create one
#			if ( $id eq "" ) {
#				$id =
#				  $parent . "-UTR"
#				  . ( $self->{genes}->{$curr_gene_name}->get_transcript($parent)
#					  ->num_UTRs() + 1 );
#			}
#			print STDERR
#			  "UTR3: $id\tTranscript: $parent\tGene: $curr_gene_name\n"
#			  if $debug;
#
#			$temp_feature =
#			  $self->{genes}->{$curr_gene_name}->get_transcript($parent)
#			  ->add_3UTR( $id, $start, $end );
#		}
#
#		
#		
#		
#
#	foreach my $ref_CDS ( @{ $self->{'exon'} } ) {
#					
#
#			my $ref_exon = $self->add_exon( $exon_id, $ref_CDS->get_start(), $ref_CDS->get_end() );
#			
#			$self->get_parent()->get_file()->add_exon_to_index( $gene_id, $self, $self->get_parent() ); 
#			
#			$ref_exon->attach_CDS($ref_CDS);
#		}
#		
#	}
#}






sub force_multi_cds {
	my $self = shift;

	# No reason to modify if a single exon transcript
	return if $self->num_exons() == 1;

	# No reason to modify if transcript already have multiple CDS, CDS > 1.
	# Or if non-coding trancript, CDS = 0
	return if ( $self->num_CDSs() != 1 );

	# No reason to modify if transcript has multple exons and only one CDS,
	# but this CDS its attached to one of the exons.
	# Meaning that this transcript has only one coding exon
	return if ( $self->{CDS}[0]->is_attached() );

	my $exon_array = $self->get_exon_array();

	# Create a copy of the single CDS
	# and delete it from Transcript
	my $single_cds = GFFCDS::copy( $self->{CDS}[0] );
	$self->delete_CDS_NO_UTR_CHANGES( $self->{CDS}[0] );

	my $single_cds_start = $single_cds->get_start();
	my $single_cds_end   = $single_cds->get_end();

	my $strand = $self->get_strand();

	foreach my $ref_exon ( @{$exon_array} ) {
		my $exon_start = $ref_exon->get_start();
		my $exon_end   = $ref_exon->get_end();

		# Skip exon if it doesn't overlap with CDS: non-coding exon
		next
		  if ( $exon_end < $single_cds_start
			|| $exon_start > $single_cds_end );

		my $new_cds_start = $exon_start;

		# Adjust new CDS start to the single CDS START if current exon
		# is the first coding exon
		if (   $single_cds_start > $exon_start
			&& $single_cds_start <= $exon_end )
		{
			$new_cds_start = $single_cds_start;
		}

		my $new_cds_end = $exon_end;

		# Adjust new CDS end to the single CDS END if current exon
		# is the LAST coding exon
		if (   $single_cds_end >= $exon_start
			&& $single_cds_end < $exon_end )
		{
			$new_cds_end = $single_cds_end;
		}

		my $temp_feature =
		  $self->add_CDS( $single_cds->get_id(), $new_cds_start, $new_cds_end,
			"." );
		$temp_feature->cpy_attrib($single_cds)

	}

	$self->attach_CDS();
	$self->sort_CDS_arr();

}

sub correct_cds_phase {
	my $self = shift;

	my $cds_sorted_array = $self->get_CDS_array_strand_based();

	my $phase_next_cds = 0;
	my $first_exon     = 1;
	foreach my $ref_cds ( @{$cds_sorted_array} ) {

		if ( $ref_cds->get_phase() ne $phase_next_cds ) {
			print STDERR "Adjusting phase on transcript "
			  . $self->{id}
			  . " CDS start: "
			  . $ref_cds->get_start()
			  . " CDS end: "
			  . $ref_cds->get_end()
			  . " current phase: "
			  . $ref_cds->get_phase()
			  . " CORRECT phase: $phase_next_cds\n";

			$ref_cds->set_phase($phase_next_cds);
		}

		# Length of the last codon encoded by this CDS feature
		my $length_last_codon = (
			abs( $ref_cds->get_start() - $ref_cds->get_end() ) + 1 -
			  $phase_next_cds ) % 3;

		$phase_next_cds = 0;
		$phase_next_cds = 1 if ( $length_last_codon == 2 );
		$phase_next_cds = 2 if ( $length_last_codon == 1 );

		$first_exon = 0;
	}

}

sub get_exon_array {
	my $self = shift;
	return $self->{exons};
}

sub get_UTR_array {
	my $self = shift;
	return $self->{UTR};
}

sub get_5prime_UTR_length {
	my $self = shift;

	my $length = 0;
	foreach my $curr_UTR ( @{ $self->{UTR} } ) {

		# type = [five_prime_utr,three_prime_utr]

		$length += $curr_UTR->get_length()
		  if $curr_UTR->get_type() eq "five_prime_utr";
	}

	return $length;
}

sub get_3prime_UTR_length {
	my $self = shift;

	my $length = 0;
	foreach my $curr_UTR ( @{ $self->{UTR} } ) {

		# type = [five_prime_utr,three_prime_utr]

		$length += $curr_UTR->get_length()
		  if $curr_UTR->get_type() eq "three_prime_utr";
	}

	return $length;
}

sub get_CDS_array {
	my $self = shift;
	return $self->{CDS};
}

sub get_CDS_length {
	my $self = shift;

	my $length = 0;
	foreach my $curr_CDS ( @{ $self->{CDS} } ) {
		$length += $curr_CDS->get_length();
	}

	return $length;
}

#sub get_CDS_length_within {
#	my $self = shift;
#
#	my ( $rng_start, $rng_end ) = @_;
#
#	my $CDS_arr  = $self->{exons};
#	my $num_CDS = scalar( @{$CDS_arr} );
#
#	return 0 if $num_CDS == 0;
#
#	my $len = 0;
#	
#	my $rngInterval = new Interval( $rng_start, $rng_end );
#
#	for ( my $CDS_ind = 0 ; $CDS_ind < ( $CDS_exons - 1 ) ; $CDS_ind++ ) {		
#		my $cdsInterval = new Interval( $CDS_arr->[$CDS_ind]->get_start(), $CDS_arr->[$CDS_ind]->get_end() );
#		
#		$len += $cdsInterval->overlap_length( $rngInterval );
#	}
#	return $len;
#}



sub get_exon_array_strand_based {
	my $self = shift;
	return $self->{exons_strand};
}

sub get_UTR_array_strand_based {
	my $self = shift;
	return $self->{UTR_strand};
}

sub get_CDS_array_strand_based {
	my $self = shift;
	return $self->{CDS_strand};
}

#sub is_identical {
#	my $self = shift;
#	my ( $other_transcript ) = @_;

#	for(my $cont_UTR = 0; $cont_UTR <= scalar( @{$self->{'UTR'}} ); $cont_UTR++ ){
#     		return 0 if not $ref_exon
#    }

#    foreach my $ref_exon ( @{$self->{exons}} ){
#      		$ref_exon->toGFF( $chrom, $self->{id}, $strand );
#    }

#    foreach my $ref_exon ( @{$self->{'CDS'}} ){
#      		$ref_exon->toGFF( $chrom, $self->{id}, $strand );
#    }
#}

sub add_exon {
	my $self = shift;
	my ( $exon_id, $start_coord, $end_coord ) = @_;
	my $temp_exon = GFFExon::new( $exon_id, $start_coord, $end_coord, $self );
	push( @{ $self->{exons} }, $temp_exon );

	my @temp = sort { $a->get_start() <=> $b->get_start() } @{ $self->{exons} };
	@{ $self->{exons} } = @temp;

	return $temp_exon;
}

sub copy_exon {
	my $self = shift;
	my ($other) = @_;

	my $temp_exon = GFFExon::copy( $other, $self );
	push( @{ $self->{exons} }, $temp_exon );

	my @temp = sort { $a->get_start() <=> $b->get_start() } @{ $self->{exons} };
	@{ $self->{exons} } = @temp;

	return $temp_exon;
}

sub delete_exon {
	my $self                 = shift;
	my ($exon_to_be_deleted) = @_;
	my $succesfuly_deleted   = 0;

	print STDERR "GFFTranscript::delete_exon ID param:" . $exon_to_be_deleted->get_id() . "\n";
	print STDERR "GFFTranscript::delete_exon Number of exons in transcript "
	  . $self->get_id() . ": "
	  . scalar( @{ $self->{exons} } ) . "\n";

	for ( my $arrInd = 0 ; $arrInd < scalar( @{ $self->{exons} } ) ; $arrInd++ )
	{
		if (
			$exon_to_be_deleted->get_id() eq $self->{exons}->[$arrInd]->get_id()
		  )
		{

			print STDERR "GFFTranscript::delete_exon ID found:"
			  . $self->{exons}->[$arrInd]->get_id() . "\n";

			# Delete upstream UTR if there is one
			my $tempUTR_UP = $self->{exons}->[$arrInd]->get_UTR("up");

			if ( $tempUTR_UP != 0 ){
				print STDERR ">>> Deleting from internal array upstream UTR "
				  . $tempUTR_UP->get_id() . "\n";
				$self->delete_UTR_from_internal_array($tempUTR_UP)
			}
			  

			# Delete downstream UTR if there is one and if it is not the same UTR as "upstream UTR"
			my $tempUTR_DOWN = $self->{exons}->[$arrInd]->get_UTR("down");

			if ( $tempUTR_DOWN != 0 && $tempUTR_DOWN != $tempUTR_UP ){
				print STDERR ">>> Deleting from internal array downstream UTR "
				  . $tempUTR_DOWN->get_id() . "\n";
				$self->delete_UTR_from_internal_array($tempUTR_DOWN)
			}
			
			# Delete CDS
			if ( $exon_to_be_deleted->has_CDS() ){
				my $tempCDS = $exon_to_be_deleted->get_CDS(); 
				print STDERR ">>> Deleting from internal array CDS "
				  . $tempCDS->get_id() . "\n";
				$self->delete_CDS_from_internal_array($tempCDS)
			}
			  

			splice @{ $self->{exons} }, $arrInd, 1;

			# Adjust transcript boundaries
			$self->adjust_boundaries_exon_based();

			# Adjust gene boundaries
			$self->get_parent()->adjust_boundaries_transcript_based();

			$succesfuly_deleted = 1;
			last;
		}
	}

	die "Unable to delete exon \'"
	  . $exon_to_be_deleted->get_id()
	  . "\' from transcript \'"
	  . $self->get_id()
	  . "\' ! Exon not found!"
	  if $succesfuly_deleted == 0;

}

sub adjust_boundaries_exon_based {
	my $self = shift;

	$self->set_start( $self->get_exon_span_start() );
	$self->set_end( $self->get_exon_span_end() );

}

sub delete_UTR_from_exon {
	my $self = shift;

	my ( $exon_with_UTR_to_be_deleted, $up_downstream ) = @_;
	my $succesfuly_deleted = 0;

	die "Unable to delete UTR from exon \'"
	  . $exon_with_UTR_to_be_deleted->get_id()
	  . "\' from transcript \'"
	  . $self->get_id()
	  . "\' ! Exon does not have an UTR!"
	  if not $exon_with_UTR_to_be_deleted->has_UTR();

	for ( my $arrInd = 0 ; $arrInd < scalar( @{ $self->{exons} } ) ; $arrInd++ )
	{
		if ( $exon_with_UTR_to_be_deleted->get_id() eq
			$self->{exons}->[$arrInd]->get_id() )
		{

			print STDERR "GFFTranscript::delete_UTR_from_exon ID param:"
			  . $exon_with_UTR_to_be_deleted->get_id() . "\n";
			print STDERR "GFFTranscript::delete_UTR_from_exon ID found:"
			  . $self->{exons}->[$arrInd]->get_id() . "\n";

			# Retrieve reference to the UTR that will be deleted
			my $tempUTR = $self->{exons}->[$arrInd]->get_UTR($up_downstream);

			die "Error. Not able to find a "
			  . $up_downstream
			  . "stream UTR on exon \'"
			  . $self->{exons}->[$arrInd]->get_id() . "\'\n"
			  if $tempUTR == 0;

			# The concept of UTR in the context of GFF is different
			# then the biological concept. The UTR is associated to an exon,
			# Example:
			# | = coding
			# : = UTR
			# - = intron
            		#
			# ||||||||||::::----------::::::::::::
			#           ^^^^ this UTR is considered internal and cannot be removed 
  
			die "Error. $up_downstream UTR "
			  . $self->{exons}->[$arrInd]->get_id() . " is internal. It cannot be removed!!'\n"
			  if $tempUTR->is_internal();
			
			# Delete UTR from Transcript object
			$self->delete_UTR($tempUTR);

			# Adjust transcript boundaries
			$self->adjust_boundaries_exon_based();

			# Adjust gene boundaries
			$self->get_parent()->adjust_boundaries_transcript_based();

			$succesfuly_deleted = 1;
			last;
		}
	}

	die "Unable to delete UTR from exon \'"
	  . $exon_with_UTR_to_be_deleted->get_id()
	  . "\' from transcript \'"
	  . $self->get_id()
	  . "\' ! Exon not found!"
	  if $succesfuly_deleted == 0;
}

sub copy_UTR {
	my $self = shift;
	my ($other) = @_;

	my $temp_utr = GFFUTR::copy( $other, $self );
	push( @{ $self->{'UTR'} }, $temp_utr );

	my @temp = sort { $a->get_start() <=> $b->get_start() } @{ $self->{'UTR'} };
	@{ $self->{'UTR'} } = @temp;

	return $temp_utr;
}

sub add_5UTR {

	my $self = shift;
	my ( $exon_id, $start_coord, $end_coord ) = @_;

	my $temp_5utr =
	  GFFUTR::new( $exon_id, $start_coord, $end_coord, "five_prime_utr",
		$self );
	push( @{ $self->{'UTR'} }, $temp_5utr );

	my @temp = sort { $a->get_start() <=> $b->get_start() } @{ $self->{'UTR'} };
	@{ $self->{'UTR'} } = @temp;

	return $temp_5utr;
}

sub delete_UTR {
	my $self                = shift;
	my ($UTR_to_be_deleted) = @_;
	my $succesfuly_deleted  = 0;

	for ( my $arrInd = 0 ; $arrInd < scalar( @{ $self->{UTR} } ) ; $arrInd++ ) {
		
		
		if ( $UTR_to_be_deleted->get_id() eq $self->{UTR}->[$arrInd]->get_id() )
		{

			my $param = $UTR_to_be_deleted->get_id();
			my $found = $self->{UTR}->[$arrInd]->get_id();
			
			print STDERR "GFFTranscript::delete_UTR ID param:" . $param . "\n";
			print STDERR "GFFTranscript::delete_UTR ID found:" . $found  . "\n";

			print STDERR "GFFTranscript::delete_UTR num UTRs:" . $self->num_UTRs() . "\n";
			

			my $parentExon = $self->{UTR}->[$arrInd]->get_parent_exon();
			$parentExon->delete_UTR( $self->{UTR}->[$arrInd] );

			print STDERR "GFFTranscript::delete_UTR num UTRs:" . $self->num_UTRs() . "\n";
			

			# Adjust transcript boundaries
			$self->adjust_boundaries_exon_based();

			# Adjust gene boundaries
			$self->get_parent()->adjust_boundaries_transcript_based();

			$succesfuly_deleted = 1;
			last;
		}
	}

	die "Unable to delete UTR \'"
	  . $UTR_to_be_deleted->get_id()
	  . "\' from transcript \'"
	  . $self->get_id()
	  . "\' ! UTR not found!"
	  if $succesfuly_deleted == 0;

}

# Delete the CDS from the exon but do not make any change
# in the UTR. !!!!! PRIVATE method !!!!!!
sub delete_CDS_NO_UTR_CHANGES {
	my $self                = shift;
	my ($CDS_to_be_deleted) = @_;
	my $succesfuly_deleted  = 0;

	for ( my $arrInd = 0 ; $arrInd < scalar( @{ $self->{CDS} } ) ; $arrInd++ ) {
		my $uu = $CDS_to_be_deleted->get_id();
		my $oo = $self->{CDS}->[$arrInd]->get_id();
		if ( $CDS_to_be_deleted->get_id() eq $self->{CDS}->[$arrInd]->get_id() )
		{

			my $parentExon = $self->{CDS}->[$arrInd]->get_parent_exon();

			# There is a possibility that the parentExon was not defined:
			# CDS is longer than exon, what prevents attachment of the CDS to
			# a parentExon. This happens in case of GFF files which the CDS
			# is a single feature.
			$parentExon->delete_CDS_NO_UTR_CHANGES( $self->{CDS}->[$arrInd] )
			  if defined $parentExon;

			splice @{ $self->{CDS} }, $arrInd, 1;

			$succesfuly_deleted = 1;
			last;
		}
	}

	die "Unable to delete CDS \'"
	  . $CDS_to_be_deleted->get_id()
	  . "\' from transcript \'"
	  . $self->get_id()
	  . "\' ! CDS not found!"
	  if $succesfuly_deleted == 0;

}

sub delete_UTR_from_internal_array {
	my $self                = shift;
	my ($UTR_to_be_deleted) = @_;
	my $succesfuly_deleted  = 0;

	for ( my $arrInd = 0 ; $arrInd < scalar( @{ $self->{UTR} } ) ; $arrInd++ ) {
		my $uu = $UTR_to_be_deleted->get_id();
		my $oo = $self->{UTR}->[$arrInd]->get_id();

		if ( $UTR_to_be_deleted->get_id() eq $self->{UTR}->[$arrInd]->get_id() )
		{

			splice @{ $self->{UTR} }, $arrInd, 1;
			$succesfuly_deleted = 1;
			last;
		}
	}

	die "Unable to delete UTR \'"
	  . $UTR_to_be_deleted->get_id()
	  . "\' from transcript \'"
	  . $self->get_id()
	  . "\' ! UTR not found!"
	  if $succesfuly_deleted == 0;

}

sub delete_CDS_from_internal_array {
	my $self                = shift;
	my ($CDS_to_be_deleted) = @_;
	my $succesfuly_deleted  = 0;

	for ( my $arrInd = 0 ; $arrInd < scalar( @{ $self->{CDS} } ) ; $arrInd++ ) {

		my $currCDS = $self->{CDS}->[$arrInd];
		
		# CDS could have the same ID. The unique identifier of a CDS is then Parent and coordinates
		if ( $CDS_to_be_deleted->get_parent()->get_id() eq $currCDS->get_parent()->get_id() &&
		     $CDS_to_be_deleted->get_start() == $currCDS->get_start() &&
		     $CDS_to_be_deleted->get_end() == $currCDS->get_end() )
		
		{
			splice @{ $self->{CDS} }, $arrInd, 1;
			$succesfuly_deleted = 1;
			last;
		}
	}

	die "Unable to delete CDS \'"
	  . $CDS_to_be_deleted->get_id()
	  . "\' from transcript \'"
	  . $self->get_id()
	  . "\' ! CDS not found!"
	  if $succesfuly_deleted == 0;

}

sub add_3UTR {

	my $self = shift;
	my ( $exon_id, $start_coord, $end_coord ) = @_;

	my $temp_3utr =
	  GFFUTR::new( $exon_id, $start_coord, $end_coord, "three_prime_utr",
		$self );

	push( @{ $self->{'UTR'} }, $temp_3utr );

	my @temp = sort { $a->get_start() <=> $b->get_start() } @{ $self->{'UTR'} };
	@{ $self->{'UTR'} } = @temp;

	return $temp_3utr;
}

sub copy_CDS {
	my $self = shift;
	my ($other) = @_;

	my $temp_cds = GFFCDS::copy( $other, $self );
	push( @{ $self->{'CDS'} }, $temp_cds );

	my @temp = sort { $a->get_start() <=> $b->get_start() } @{ $self->{'CDS'} };
	@{ $self->{'CDS'} } = @temp;

	return $temp_cds;
}

sub add_CDS {

	my $self = shift;
	my ( $exon_id, $start_coord, $end_coord, $phase ) = @_;

	my $temp_cds =
	  GFFCDS::new( $exon_id, $start_coord, $end_coord, $phase, $self );
	push( @{ $self->{'CDS'} }, $temp_cds );

	my @temp = sort { $a->get_start() <=> $b->get_start() } @{ $self->{'CDS'} };
	@{ $self->{'CDS'} } = @temp;

	return $temp_cds;
}

sub get_exon_span_start {
	my $self = shift;
	return $self->{exons}[0]->get_start();
}

sub get_exon_span_end {
	my $self = shift;
	return $self->{exons}[ $self->num_exons() - 1 ]->get_end();
}

sub set_id {
	my $self = shift;

	my ($new_id) = @_;

	$self->{id} = $new_id;
}

sub get_id {
	my $self = shift;
	return $self->{id};
}

sub get_start {
	my $self = shift;
	return $self->{start};
}

sub get_end {
	my $self = shift;
	return $self->{end};
}

sub get_length {
	my $self = shift;
	return $self->{end} - $self->{start};
}

sub set_start {
	my $self = shift;
	my ($start_coord) = @_;
	$self->{start} = $start_coord;
}

sub set_end {
	my $self = shift;
	my ($end_coord) = @_;
	return $self->{end} = $end_coord;
}

sub num_exons {
	my $self = shift;
	if ( defined( $self->{exons} ) ) {
		return scalar( @{ $self->{exons} } );
	}
	else {
		return 0;
	}
}

sub num_CDSs {
	my $self = shift;
	if ( defined( $self->{CDS} ) ) {
		return scalar( @{ $self->{CDS} } );
	}
	else {
		return 0;
	}
}

sub num_UTRs {
	my $self = shift;
	if ( defined( $self->{UTR} ) ) {
		return scalar( @{ $self->{UTR} } );
	}
	else {
		return 0;
	}
}

sub get_first_CDS_strand_based {
	my $self = shift;
	if ( defined( $self->{CDS_strand} ) ) {
		return $self->{CDS_strand}[0];
	}
	else {
		return 0;
	}
}

sub get_last_CDS_strand_based {
	my $self = shift;
	if ( defined( $self->{CDS_strand} ) ) {
		return $self->{CDS_strand}[ scalar( @{ $self->{CDS_strand} } ) - 1 ];
	}
	else {
		return 0;
	}
}


sub get_first_exon {
	my $self = shift;
	if ( defined( $self->{exons} ) ) {
		return $self->{exons}[0];
	}
	else {
		return 0;
	}
}

sub get_first_exon_strand_based {
	my $self = shift;
	if ( defined( $self->{exons_strand} ) ) {
		return $self->{exons_strand}[0];
	}
	else {
		return 0;
	}
}

sub get_last_exon {
	my $self = shift;
	if ( defined( $self->{exons} ) ) {
		return $self->{exons}[ scalar( @{ $self->{exons} } ) - 1 ];
	}
	else {
		return 0;
	}
}

sub get_last_exon_strand_based {
	my $self = shift;
	if ( defined( $self->{exons_strand} ) ) {
		return $self->{exons_strand}[ scalar( @{ $self->{exons} } ) - 1 ];
	}
	else {
		return 0;
	}
}

sub get_first_UTR {
	my $self = shift;
	if ( defined( $self->{UTR} ) ) {
		return ${ $self->{UTR} }[0];
	}
	else {
		return 0;
	}
}

sub get_last_UTR {
	my $self = shift;
	if ( defined( $self->{UTR} ) ) {
		return $self->{UTR}[ scalar( @{ $self->{UTR} } ) - 1 ];
	}
	else {
		return 0;
	}
}

sub get_exon_by_name {
	my $self = shift;

	my ($exon_id) = @_;

	return $self->get_exon($exon_id);
}

#######################
# Poor implementation
#
# Transcript should store exons in a array and also in a hash
# then the loop will be not necessary

sub get_exon {
	my $self = shift;

	my ($exon_id) = @_;

	foreach my $ref_exon ( @{ $self->{exons} } ) {
		return $ref_exon if ( $ref_exon->get_id() eq $exon_id );
	}

	die "GFFTransript::get_exon_by_name(): There are no exons with id ="
	  . $exon_id . "\n";
}

sub get_exon_by_coord {
	my $self = shift;

	my ($coord) = @_;

	foreach my $ref_exon ( @{ $self->{exons} } ) {
		return $ref_exon if ( $ref_exon->get_start() <= $coord && $ref_exon->get_end() >= $coord );
	}

	return 0;
}


sub get_exon_by_start_end_coord {
	my $self = shift;

	my ($start, $end) = @_;

	foreach my $ref_exon ( @{ $self->{exons} } ) {
		return $ref_exon if ( $ref_exon->get_start() == $start && $ref_exon->get_end() == $end );
	}

	return 0;
}


sub get_exon_with_same_start {
	my $self = shift;

	my ($other_exon_ref) = @_;

	foreach my $ref_exon ( @{ $self->{exons} } ) {
		return $ref_exon if ( $ref_exon->has_same_start($other_exon_ref) );
	}

	return 0;
}

sub get_exon_with_same_end {
	my $self = shift;

	my ($other_exon_ref) = @_;

	foreach my $ref_exon ( @{ $self->{exons} } ) {
		return $ref_exon if ( $ref_exon->has_same_end($other_exon_ref) );
	}

	return 0;
}

sub get_compatible_exon {
	my $self = shift;

	my ($other_exon_ref) = @_;

	foreach my $ref_exon ( @{ $self->{exons} } ) {
		return $ref_exon if ( $ref_exon->is_compatible_to($other_exon_ref) );
	}

	return 0;
}

# feature_type = [exon,five_prime_utr,three_prime_utr,CDS]
sub num_child_features {
	my $self = shift;

	my ($feature_type) = @_;

	if ( defined( $self->{$feature_type} ) ) {
		return scalar( @{ $self->{$feature_type} } );
	}
	else {
		return 0;
	}
}

sub push {
	my $self = shift;

	my ( $other, $start, $end ) = @_;

	foreach my $cexon ( @{ $other->{exons} } ) {
		if ( $cexon->get_start() >= $start && $cexon->get_end() <= $end ) {
			my $tmp_exon = $self->copy_exon($cexon);

			if ( $cexon->has_UTR() ) {

				# Up
				my $up_utr = $cexon->get_UTR("up");

				if ( $up_utr != 0 ) {
					my $tmp_utr = $self->copy_UTR($up_utr);
					$tmp_exon->attach_UTR($tmp_utr);
				}

				# Down
				my $down_utr = $cexon->get_UTR("down");

				if ( $down_utr != 0 ) {
					my $tmp_utr = $self->copy_UTR($down_utr);
					$tmp_exon->attach_UTR($tmp_utr);
				}

			}

			if ( $cexon->has_CDS() ) {
				my $tmp_cds = $self->copy_CDS( $cexon->get_CDS() );
				$tmp_exon->attach_CDS($tmp_cds);
			}

		}
	}
	
	$self->adjust_boundaries_exon_based();
	$self->get_parent()->adjust_boundaries_transcript_based();

	return;
}

sub toGFF {
	my $self = shift;

	my ( $chrom, $parent, $strand ) = @_;

	my $str =
	    $chrom . "\t.\t"
	  . $self->{transcript_type} . "\t"
	  . $self->{start} . "\t"
	  . $self->{end} . "\t.\t"
	  . $strand . "\t.\t" . "ID="
	  . $self->{id} . ";"
	  . "Parent=$parent;\n";

	foreach my $ref_exon ( @{ $self->{'UTR'} } ) {
		$str .= $ref_exon->toGFF( $chrom, $self->{id}, $strand );
	}

	foreach my $ref_exon ( @{ $self->{exons} } ) {
		$str .= $ref_exon->toGFF( $chrom, $self->{id}, $strand );
	}

	foreach my $ref_exon ( @{ $self->{'CDS'} } ) {
		$str .= $ref_exon->toGFF( $chrom, $self->{id}, $strand );
	}

	return $str;
}

sub toGTF {
	my $self = shift;

	my ( $chrom, $parent, $strand ) = @_;
		
	my $str;
	foreach my $ref_exon ( @{ $self->{'UTR'} } ) {
		#print STDERR "IN UTR\n";
		#getc();
		$str .= $ref_exon->toGTF( $chrom, $self->{id}, $strand );
	}

	foreach my $ref_exon ( @{ $self->{'CDS'} } ) {
		#print STDERR "IN CDS\n";
		#getc();
		$str .= $ref_exon->toGTF( $chrom, $self->{id}, $strand );
	}

	my $transcript_id = $self->get_id();
	my $gene_id = $self->get_parent()->get_id();
	
	$str =~ s/\n/\tgene_id \"$gene_id\"; transcript_id \"$transcript_id\";\n/g;

	return $str;
}


sub set_attribute {
	my $self = shift;
	my ( $key, $value ) = @_;

	return $self->{attrib}{$key} = $value;
}

sub get_attribute {
	my $self = shift;
	my ($key) = @_;
	return $self->{attrib}{$key};
}

sub get_strand {
	my $self = shift;
	return $self->{parent}->get_strand();
}

sub get_chrom {
	my $self = shift;
	return $self->{parent}->get_chrom();
}

return 1;
