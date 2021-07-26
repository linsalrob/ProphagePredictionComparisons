package GFFCDS;
use strict 'vars';
use strict 'refs';

require GFFTranscript;

# use overload '==' => \&compare;

sub new {
	my ( $exon_name, $start_coord, $end_coord, $phase, $parent ) = @_;
	my $self = {
		id     => $exon_name,
		start  => $start_coord,
		end    => $end_coord,
		phase  => $phase,
		parent => $parent
	};

	bless $self, GFFCDS;

	return $self;
}

sub copy {
	my ( $ref, $parent ) = @_;

	my $self;

	if ( defined $parent ) {
		$self = {
			id     => $ref->{id},
			start  => $ref->{start},
			end    => $ref->{end},
			phase  => $ref->{phase},
			parent => $parent
		};
	}
	else {
		$self = {
			id     => $ref->{id},
			start  => $ref->{start},
			end    => $ref->{end},
			phase  => $ref->{phase},
			parent => $ref->{parent}
		};
	}

	bless $self, GFFCDS;

	# Copying other atrributes
	$self->cpy_attrib($ref);
	return $self;
}

sub cpy_attrib {
	my $self = shift;
	my ($ref) = @_;

	foreach my $ckey ( keys %{ $ref->{attrib} } ) {
		$self->{attrib}{$ckey} = $ref->{attrib}{$ckey};
	}
}

sub get_parent {
	my $self = shift;
	return $self->{parent};
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
	return $self->{end} - $self->{start} + 1;
}

sub get_phase {
	my $self = shift;
	return $self->{phase};
}

sub set_phase {
	my $self = shift;
	my ($phase) = @_;
	$self->{phase} = $phase;
}

sub set_start {
	my $self = shift;
	my ($coord_start) = @_;

	$self->{start} = $coord_start;
}

sub set_end {
	my $self = shift;
	my ($coord_end) = @_;

	$self->{end} = $coord_end;
}

sub set_parent_exon {
	my $self = shift;
	my ($ref_exon) = @_;

	$self->{exon} = $ref_exon;
}

sub get_parent_exon {
	my $self = shift;
	return $self->{exon};
}

sub is_attached {
	my $self = shift;
	return defined $self->{exon};
}

sub toGFF {
	my $self = shift;

	my ( $chrom, $parent, $strand ) = @_;

	my $str =
	    $chrom . "\t.\t" . "CDS" . "\t"
	  . $self->{start} . "\t"
	  . $self->{end} . "\t.\t"
	  . $strand . "\t"
	  . $self->{phase} . "\t" . "ID="
	  . $self->{id} . ";"
	  . "Parent=$parent;\n";
	return $str;
}

sub toGTF {
	my $self = shift;

	my ( $chrom, $parent, $strand ) = @_;

	my $str;
	
	#print STDERR "INSIDE CDS\n";
	#getc();

	# If the CDS is the first one in the transcript 
	if( $self->is_first_strand_based() ){
		my $start_start_codon;
		my $end_start_codon;
		
	#print STDERR "FIRST CDS\n";
	#getc();
		

		if ( $strand eq '+' ) {
			$start_start_codon = $self->{start};
			$end_start_codon =  $self->{start} + 2;
		}
		else {
			$start_start_codon = $self->{end} - 2;
			$end_start_codon = $self->{end};
		}
		
		$str .=
		    $chrom . "\t.\t" . "start_codon" . "\t"
		  . $start_start_codon . "\t"
		  . $end_start_codon . "\t.\t"
		  . $strand . "\t"
		  . $self->{phase} . "\n";



	# If the CDS is the last one in the transcript the STOP codon
	# should not be part of it
	}
	
	if ( $self->is_last_strand_based() ) {

	#print STDERR "LAST CDS\n";
	#getc();
		

		my $start = $self->{start};
		my $end   = $self->{end};

		if ( $strand eq '+' ) {
			$end -= 3;
		}
		else {
			$start += 3;
		}

		if( $start < $end ){
		$str .=
		    $chrom . "\t.\t" . "CDS" . "\t" 
		  . $start . "\t" 
		  . $end . "\t.\t"
		  . $strand . "\t"
		  . $self->{phase} . "\n";
		}
		  
		  
		my $start_stop_codon;
		my $end_stop_codon;

		if ( $strand eq '+' ) {
			$start_stop_codon = $self->{end} - 2;
			$end_stop_codon =  $self->{end};
		}
		else {
			$start_stop_codon = $self->{start};
			$end_stop_codon = $self->{start} + 2;
		}
		
		$str .=
		    $chrom . "\t.\t" . "stop_codon" . "\t"
		  . $start_stop_codon . "\t"
		  . $end_stop_codon . "\t.\t"
		  . $strand . "\t"
		  . $self->{phase} . "\n";
		  

	} else {
		
	#print STDERR "REGULAR CDS\n";
	#getc();
		
		$str .=
		    $chrom . "\t.\t" . "CDS" . "\t"
		  . $self->{start} . "\t"
		  . $self->{end} . "\t.\t"
		  . $strand . "\t"
		  . $self->{phase} . "\n";

	}

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

sub is_last_strand_based {
	my $self = shift;

	return $self eq $self->{parent}->get_last_CDS_strand_based();
}

sub is_first_strand_based {
	my $self = shift;

	return $self eq $self->{parent}->get_first_CDS_strand_based();
}

return 1;
