package GFFExon;
use strict 'vars';
use strict 'refs';

require GFFTranscript;

use GFFCDS;
use GFFUTR;

my $debug = 1;

# use overload '==' => \&compare;

sub new {
	my ( $exon_name, $start_coord, $end_coord, $parent ) = @_;
	my $self = { id => $exon_name,
		     start => $start_coord,
		     end => $end_coord,
		     parent => $parent };
			  		     
	bless $self, GFFExon;
	
	return $self;
}

sub copy {
	my ( $ref, $parent ) = @_;


	my $self = {
			id     => $ref->{id},
			start  => $ref->{start},
			end    => $ref->{end},
			parent => $parent
		};

	bless $self, GFFExon;

	# Copying other atrributes
	$self->cpy_attrib($ref);
	return $self;
}

sub cpy_attrib{
	my $self = shift;
	my ( $ref ) = @_;
	
	foreach my $ckey ( keys %{$ref->{attrib}} ){
		$self->{attrib}{$ckey} = $ref->{attrib}{$ckey}
	}
}


sub overlaps {
	my $self = shift;
	my ($other) = @_;
	
	return 0 if $self->get_chrom() ne $other->get_chrom();
	return 0 if $self->get_strand() ne $other->get_strand();
	
	if( ( $self->get_start() >= $other->get_start && $self->get_start() <= $other->get_end() ) ||
	    ( $other->get_start() >= $self->get_start && $other->get_start() <= $self->get_end() ) ){
	    	return 1;
	    }else{
	    	return 0;
	    }
}

sub overlaps_op_strand {
	my $self = shift;
	my ($other) = @_;
	
	return 0 if $self->get_chrom() ne $other->get_chrom();
	return 0 if $self->get_strand() eq $other->get_strand();
	
	if( ( $self->get_start() >= $other->get_start && $self->get_start() <= $other->get_end() ) ||
	    ( $other->get_start() >= $self->get_start && $other->get_start() <= $self->get_end() ) ){
	    	return 1;
	    }else{
	    	return 0;
	    }
}



sub identical {
	my $self = shift;
	my ($other) = @_;
	
	return 0 if $self->get_chrom() ne $other->get_chrom();
	return 0 if $self->get_strand() ne $other->get_strand();
	
	if( $self->get_start() == $other->get_start && $self->get_end() == $other->get_end() ){
	    	return 1;
	    }else{
	    	return 0;
	    }
}

sub get_strand{
	my $self = shift;
	return $self->{parent}->get_strand();	
}

sub get_chrom{
	my $self = shift;
	return $self->{parent}->get_chrom();	
}


sub get_parent {
	my $self = shift;
	return $self->{parent};
}

sub is_last{
	my $self = shift;
	
	#print "GFFExonUTRCDS::is_last   Parent ID " . $self->{parent}->get_id() . "\n";
	#print "Last exon:" . $self->{parent}->get_last_exon()->get_id() . "  Current exon:" . $self->{id} . "\n";
	
	return $self->{id} eq $self->{parent}->get_last_exon()->get_id();
}

sub is_last_exon_strand_based{
	my $self = shift;
	
	return $self->{id} eq $self->{parent}->get_last_exon_strand_based()->get_id();
}	

sub is_first_exon_strand_based{
	my $self = shift;
	
	return $self->{id} eq $self->{parent}->get_first_exon_strand_based()->get_id();
}	
	

sub is_first{
	my $self = shift;
	
	#print "GFFExonUTRCDS::is_first   Parent ID " . $self->{parent}->get_id() . "\n";
	#print "First exon:" . $self->{parent}->get_first_exon()->get_id() . "  Current exon:" . $self->{id} . "\n";
	
	return $self->{id} eq $self->{parent}->get_first_exon()->get_id();
}


sub has_same_start{
	my $self = shift;

	my ($other_exon_ref) = @_;

	return 1 if( ( $self->{start} == $other_exon_ref->get_start() )  );
		
	return 0;
}

sub has_same_end{
	my $self = shift;

	my ($other_exon_ref) = @_;

	return 1 if( ( $self->{end} == $other_exon_ref->get_end() )  );
		
	return 0;
}


sub is_compatible_to{
	my $self = shift;

	my ($other_exon_ref) = @_;

	# Exons are compatible if:	
	# - They are the last exons and start coordinates are identical
	# - They are the first exons and end coordinates are identical
	# - If they are neither the first nor the last exons but the coordinates are identical
	# - Both exons are the single exon in the gene and they have identical CDS

	# For the time being we are going compare only the same exon across different GFF files
	# Therefore, if the exon id is different the exons are automatically incompatible	



	# If both exons belong to a transcript with only one exon
	# and the CDS are identical, those exons are considered 
	# compatible
	if ( $other_exon_ref->get_parent()->num_exons() == 1 &&
	               $self->get_parent()->num_exons() == 1 ){
	  my $other_cds = $other_exon_ref->get_CDS();
	  my $self_cds  = $other_exon_ref->get_CDS();
	  
	  if( $other_cds->get_start() == $other_cds->get_start() &&
	      $other_cds->get_end()   == $other_cds->get_end() ){
	      return 1
	  }else{
	  	  return 0;
	  }	  
	}
	
	return 1 if( ( $self->{start} == $other_exon_ref->get_start() ) &&
		     ( $self->{end} == $other_exon_ref->get_end() )  );
	
	return 1 if( $self->is_last() && $other_exon_ref->is_last() &&
		     ( $self->{start} == $other_exon_ref->get_start() ) );

	return 1 if( $self->is_first() && $other_exon_ref->is_first() &&
		     ( $self->{end} == $other_exon_ref->get_end() ) );

	
	return 0;
}

sub get_id {
	my $self = shift;
	return $self->{id};
}

sub set_id {
	my $self = shift;
	
	my ($new_id) = @_;
	
	$self->{id} = $new_id;
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
	my ( $coord_start ) = @_;
	
	$self->{start} = $coord_start;
	
	foreach my $currUTR ( @{$self->{UTR}} ){  
		$currUTR->set_start( $coord_start ) if( $currUTR->get_start() < $coord_start );
		if( $currUTR->get_length() <= 0 ){
			$self->get_parent()->delete_UTR_from_internal_array( $currUTR );
			$self->dettach_UTR( $currUTR );
		}			 
	}
	
	if( defined $self->{CDS} ){	  
		$self->{CDS}->set_start( $coord_start) if( $self->{CDS}->get_start() < $coord_start );
	}  
}
	
sub set_end {
	my $self = shift;
	my ( $coord_end ) = @_;
	
	$self->{end} = $coord_end;
		
	foreach my $currUTR ( @{$self->{UTR}} ){  
		$currUTR->set_end( $coord_end) if( $currUTR->get_end() > $coord_end );		
		if( $currUTR->get_length() <= 0 ){
			$self->get_parent()->delete_UTR_from_internal_array( $currUTR );
			$self->dettach_UTR( $currUTR );
		}			 		
	}
	
	if( defined $self->{CDS} ){
		$self->{CDS}->set_end( $coord_end ) if( $self->{CDS}->get_end() > $coord_end );
	}  	
}

sub toGFF {
	my $self = shift;
	
	my ($chrom, $parent, $strand ) = @_;
	
	my $str = $chrom . "\t.\t" . "exon" . "\t" . 
	      $self->{start} . "\t" . 
	      $self->{end} . "\t.\t" . 
	      $strand . "\t.\t" . 
	      "ID=" . $self->{id} . ";" .
      	      "Parent=$parent;\n";
      	      
    return $str;	
}

sub attach_UTR{
	my $self = shift;

	my ($utr_ref) = @_;
	
	push( @{$self->{UTR}}, $utr_ref );

	$utr_ref->set_parent_exon( $self );
}

sub dettach_UTR{
	my $self = shift;

	my ($param_UTR) = @_;
	
	my $succesfuly_dettached = 0;
	
	for( my $arrInd = 0; $arrInd < scalar( @{$self->{UTR}} ); $arrInd++ ){
		if( $param_UTR->get_id() eq $self->{UTR}->[ $arrInd ]->get_id() ){
			
			print STDERR "GFFExon::dettach_UTR ID param:" . $param_UTR->get_id() . "\n";
			print STDERR "GFFExon::dettach_UTR ID found:" . $self->{UTR}->[ $arrInd ]->get_id()  . "\n";
						
			die "Something wrong! Exon " . $self->get_id() . " was reported as coding" .
			    " exon, but there are no CDS attached!"  if( not defined( $self->{CDS} ) ); 
									
			splice @{$self->{UTR}}, $arrInd, 1;			
			
			$succesfuly_dettached = 1;
			last;
		}	
	}	
	
	die "Unable to delete UTR from exon \'" . $self->get_id() . 
	    "\' ! UTR not found!" 
	    if $succesfuly_dettached == 0;  
	
}

# Delete the CDS from the exon but do not make any change
# in the UTR. !!!!! PRIVATE method !!!!!!
sub delete_CDS_NO_UTR_CHANGES{
	my $self = shift;
	my ($param_CDS) = @_;

	my $succesfuly_deleted = 0;
		
	die "Unable to delete CDS from exon \'" . $self->get_id() . 
	    "\'! Exon does not have an CDS!" 
	    if not $self->has_CDS();  

	die "Unable to delete CDS from exon \'" . $self->get_id() . 
	    "\'! Exon is non-coding!" 
	    if $self->is_non_coding();  

	print STDERR "Removing CDS: " . $self->get_id() . "\n" if $debug;
		
	for( my $arrInd = 0; $arrInd < scalar( @{$self->{CDS}} ); $arrInd++ ){
		if( $param_CDS->get_id() eq $self->{CDS}->[ $arrInd ]->get_id() ){
			
			print STDERR "GFFExon::delete_CDS_NO_UTR_CHANGES ID param:" . $param_CDS->get_id() . "\n";
			print STDERR "GFFExon::delete_CDS_NO_UTR_CHANGES ID found:" . $self->{CDS}->[ $arrInd ]->get_id()  . "\n";
						
			die "Something wrong! Exon " . $self->get_id() . " was reported as coding" .
			    " exon, but there are no CDS attached!"  if( not defined( $self->{CDS} ) ); 
									
			splice @{$self->{CDS}}, $arrInd, 1;			
			
			$succesfuly_deleted = 1;
			last;
		}	
	}

	die "Unable to delete CDS from exon \'" . $self->get_id() . 
	    "\' ! CDS not found!" 
	    if $succesfuly_deleted == 0;  

}


sub delete_UTR{
	my $self = shift;
	my ($param_UTR) = @_;

	my $succesfuly_deleted = 0;
		
	die "Unable to delete UTR from exon \'" . $self->get_id() . 
	    "\'! Exon does not have an UTR!" 
	    if not $self->has_UTR();  

	die "Unable to delete UTR from exon \'" . $self->get_id() . 
	    "\'! Exon is non-coding!" 
	    if $self->is_non_coding();  

		
	for( my $arrInd = 0; $arrInd < scalar( @{$self->{UTR}} ); $arrInd++ ){
		if( $param_UTR->get_id() eq $self->{UTR}->[ $arrInd ]->get_id() ){
			
			print STDERR "GFFExon::delete_UTR ID param:" . $param_UTR->get_id() . "\n";
			print STDERR "GFFExon::delete_UTR ID found:" . $self->{UTR}->[ $arrInd ]->get_id()  . "\n";
						
			die "Something wrong! Exon " . $self->get_id() . " was reported as coding" .
			    " exon, but there are no CDS attached!"  if( not defined( $self->{CDS} ) ); 
						
			# If UTR on upstream region of the exon
			if ( $self->{UTR}->[ $arrInd ]->get_start() == $self->get_start() ){
				$self->set_start( $self->{CDS}->get_start );
			}elsif( $self->{UTR}->[ $arrInd ]->get_end() == $self->get_end() ){
				$self->set_end( $self->{CDS}->get_end );
			}else{
				die "Something wrong on exon " . $self->get_id() . " ! UTR seems to be in the core of the EXON";
			} 

			print STDERR "GFFExon::delete_UTR Transcript numUTRs:" . $self->get_parent()->num_UTRs . "\n";
			
			
			$succesfuly_deleted = 1;
			last;
		}	
	}

	die "Unable to delete UTR from exon \'" . $self->get_id() . 
	    "\' ! UTR not found!" 
	    if $succesfuly_deleted == 0;  

}


sub is_non_coding {
	my $self = shift;
	my ( $exon_param ) = @_;
	
	return 0 if not $self->has_UTR();
	
	foreach my $currUTR ( @{$self->{UTR}} ){  

		return 1 if( $currUTR->get_start == $self->get_start() &&
				     $currUTR->get_end   == $self->get_end()  );
		
	}
	
	return 0;
}

sub has_UTR{
	my $self = shift;

	return 1 if ( defined( $self->{UTR} ) && scalar( @{$self->{UTR}} ) != 0 );
	return 0;
}	

sub has_up_UTR{
	my $self = shift;		
	foreach my $currUTR ( @{$self->{UTR}} ){
		
		    print STDERR "GFFExon::has_down_utr curr. UTR start: " . $currUTR->get_start() . "\n";
		    print STDERR "GFFExon::has_down_utr curr. EXON start: " . $self->get_start() . "\n";		    
		  		
			return 1 if $currUTR->get_start() == $self->get_start();
	}	
	return 0;
}	

sub has_down_UTR{
	my $self = shift;
	foreach my $currUTR ( @{$self->{UTR}} ){  	
		    print STDERR "GFFExon::has_down_utr curr. UTR end: " . $currUTR->get_end() . "\n";
		    print STDERR "GFFExon::has_down_utr curr. EXON end: " . $self->get_end() . "\n";		    
			return 1 if $currUTR->get_end() == $self->get_end();
	}	
	return 0;
}	


sub attach_CDS{
	my $self = shift;

	my ($cds_ref) = @_;	
	$self->{CDS} = $cds_ref;
	$cds_ref->set_parent_exon( $self );
}

sub get_UTR{
	my $self = shift;
	
	my ($up_downstream) = @_;
	
	foreach my $currUTR ( @{$self->{UTR}} ){  		
		if( $up_downstream eq "up" ){
			return $currUTR if $currUTR->get_start() == $self->get_start()
		}elsif( $up_downstream eq "down" ){
			return $currUTR if $currUTR->get_end() == $self->get_end()
		}
	}	
	
	return 0;
}

sub get_CDS{
	my $self = shift;	
	return $self->{CDS};
}

sub has_CDS{
	my $self = shift;
	return 1 if defined( $self->{CDS} );
	return 0;
}	



# sub compare{
	# my $self = shift;
	# 
	# my ( $other ) = @_;
	# 
	# return 0 if( 	$self->{start} == $other->{start} && 
			# $self->{end} == $other->{end}	  && 
			# $self->{chrom} == $other->{end}	  && 
	# 
# }

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

	
return 1;
