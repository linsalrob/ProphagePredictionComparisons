package GFFUTR;
use strict 'vars';
use strict 'refs';

require GFFTranscript;
require GFFExon;

# use overload '==' => \&compare;
# type = [five_prime_utr,three_prime_utr]	
sub new {
	my ( $exon_name, $start_coord, $end_coord, $type, $parent ) = @_;
	my $self = { id => $exon_name,
		     	 start => $start_coord,
		     	 end => $end_coord,
		     	 type => $type,		     
		     	 parent => $parent };
			  
	bless $self, GFFUTR;
	
	return $self;
}

sub copy {
	my ( $ref, $parent ) = @_;


	my $self = {
			id     => $ref->{id},
			start  => $ref->{start},
			end    => $ref->{end},
			type   => $ref->{type},
			parent => $parent
		};

	bless $self, GFFUTR;

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

			# The concept of UTR in the context of GFF is different
			# then the biological concept. The UTR is associated to an exon,
			# Example:
			# | = coding
			# : = UTR
			# - = intron
            		#
			# ||||||||||::::----------::::::::::::
			#           ^^^^ this UTR is considered internal and cannot be removed 


sub is_internal {
    my $self = shift;

    my $transcript = $self->{parent}->get_parent();
	return 0 if ( $self->get_start() == $transcript->get_start() || $self->get_end() == $transcript->get_end() );
	return 1;
}

# type = [five_prime_utr,three_prime_utr]	
sub get_type {
	my $self = shift;
	return $self->{type};
}

sub get_type_GTF {
	my $self = shift;
	return '5UTR' if $self->{type} eq 'five_prime_utr';
	return '3UTR' if $self->{type} eq 'three_prime_utr';
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

sub set_parent_exon{
	my $self = shift;
	my ( $ref_exon ) = @_;
	
	$self->{exon} = $ref_exon; 
}

sub get_parent_exon{
	my $self = shift;
	return $self->{exon}; 
}


sub toGFF {
	my $self = shift;
	
	my ($chrom, $parent, $strand ) = @_;
	
	my $str = $chrom . "\t.\t" . $self->{type} . "\t" . 
	      $self->{start} . "\t" . 
	      $self->{end} . "\t.\t" . 
	      $strand . "\t.\t" . 
	      "ID=" . $self->{id} . ";" .
      	      "Parent=$parent;\n";
      	      
    return $str;	
}

sub toGTF {
	my $self = shift;
	
	#print STDERR "UTR\n";
	#getc();
	
	my ($chrom, $parent, $strand ) = @_;
	
	my $str = $chrom . "\t.\t" . $self->get_type_GTF() . "\t" . 
	      $self->{start} . "\t" . 
	      $self->{end} . "\t.\t" . 
	      $strand . "\t.\n";
      	      
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
	
return 1;
