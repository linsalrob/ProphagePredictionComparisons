package GFFExonUTRCDS;
use strict 'vars';
use strict 'refs';

require GFFLib::GFFTranscript;

# use overload '==' => \&compare;

sub new {
# type = [exon,five_prime_utr,three_prime_utr,CDS]	
	my ( $exon_name, $start_coord, $end_coord, $type, $parent ) = @_;
	my $self = { id => $exon_name,
		     start => $start_coord,
		     end => $end_coord,
		     type => $type,
		     parent => $parent };
			  
	bless $self, GFFExonUTRCDS;
	
	return $self;
}

sub get_parent {
	my $self = shift;
	return $self->{parent};
}

# type = [exon,five_prime_utr,three_prime_utr,CDS]	
sub is_last{
	my $self = shift;
	
	print "GFFExonUTRCDS::is_last   Parent ID " . $self->{parent}->get_id() . "\n";
	print "Last exon:" . $self->{parent}->get_last( $self->{type} )->get_id() . "  Current:" . $self->{id} . "\n";
	
	return $self->{id} eq $self->{parent}->get_last( $self->{type} )->get_id();
}

sub is_first{
	my $self = shift;
	
	print "GFFExonUTRCDS::is_first   Parent ID " . $self->{parent}->get_id() . "\n";
	print "First exon:" . $self->{parent}->get_first( $self->{type} )->get_id() . "  Current:" . $self->{id} . "\n";
	
	return $self->{id} eq $self->{parent}->get_first( $self->{type} )->get_id();
}

sub is_compatible_to{
	my $self = shift;

	my ($other_exon_ref) = @_;

	# Exons are compatible if:	
	# - They are the last exons and start coordinates are identical
	# - They are the first exons and end coordinates are identical
	# - If they are neither the first nor the last exons but the coordinates are identical

	# For the time being we are going compare only the same exon across different GFF files
	# Therefore, if the exon id is different the exons are automatically incompatible	
	
	return 1 if( ( $self->{start} == $other_exon_ref->get_start() ) &&
		     ( $self->{end} == $other_exon_ref->get_end() )  );
	
	return 1 if( $self->is_last_exon() && $other_exon_ref->is_last_exon() &&
		     ( $self->{start} == $other_exon_ref->get_start() ) );

	return 1 if( $self->is_first_exon() && $other_exon_ref->is_first_exon() &&
		     ( $self->{end} == $other_exon_ref->get_end() ) );

	
	return 0;
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

sub set_start {
	my $self = shift;
	my ( $coord_start ) = @_;
	
	return $self->{start} = $coord_start;
}
	
sub set_end {
	my $self = shift;
	my ( $coord_end ) = @_;
	
	return $self->{end} = $coord_end;
}

sub toGFF {
	my $self = shift;
	
	my ($chrom, $parent, $strand ) = @_;
	
	print $chrom . "\t.\t" . $self->{type} . "\t" . 
	      $self->{start} . "\t" . 
	      $self->{end} . "\t.\t" . 
	      $strand . "\t.\t" . 
	      "ID=" . $self->{id} . "; " .
      	      "Parent=$parent;\n";	
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
