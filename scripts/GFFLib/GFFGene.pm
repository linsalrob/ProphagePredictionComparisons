package GFFGene;
use strict 'vars';
use strict 'refs';

use GFFTranscript;

sub new {
	my ( $gene_id, $gene_name, $chrom, $start_coord, $end_coord, $strand, $file ) = @_;
	my $self = {
		id     => $gene_id,
		name   => $gene_name,
		chrom  => $chrom,
		start  => $start_coord,
		end    => $end_coord,
		strand => $strand,
		file   => $file
	};
	bless $self, GFFGene;

	return $self;
}

sub get_file{
	my $self = shift;
	return $self->{file};
}


sub set_filename {
	my $self = shift;
	my ($filename) = @_;
	$self->{filename} = $filename;
}

sub get_filename {
	my $self = shift;
	return $self->{filename};
}

sub overlaps {
	my $self = shift;

	# $strand_based = [op_strand, same_strand, strandness]
	my ( $other, $strand_based, $buffer ) = @_;

	return 0 if $self->get_chrom() ne $other->get_chrom();
	return 0 if $self->get_strand() ne $other->get_strand() && $strand_based eq "same_strand";
	return 0 if $self->get_strand() eq $other->get_strand() && $strand_based eq "op_strand";

	if (
		(
			   $self->get_start() >= ( $other->get_start - $buffer )
			&& $self->get_start() <= ( $other->get_end() + $buffer )
		)
		|| (   $other->get_start() >= ( $self->get_start - $buffer )
			&& $other->get_start() <= ( $self->get_end() + $buffer ) )
	  )
	{
		return 1;
	}
	else {
		return 0;
	}
}

sub add_transcript {
	my $self = shift;
	my ( $transcript_name, $transcript_type, $start_coord, $end_coord ) = @_;

	$self->{transcript}->{$transcript_name} =
	  GFFTranscript::new( $transcript_name, $transcript_type, $start_coord, $end_coord, $self );

	return $self->{transcript}->{$transcript_name};
}

sub delete_transcript {
	my $self = shift;
	my ($id) = @_;
	undef $self->{transcript}->{$id};
	delete $self->{transcript}->{$id};
	
	$self->adjust_boundaries_transcript_based();
}

sub get_name {
	my $self = shift;
	return $self->{name};
}

sub set_name {
	my $self = shift;

	my ($new_name) = @_;
	$self->{name} = $new_name;
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

sub set_description {
	my $self = shift;
	my ($description) = @_;
	return $self->{description} = $description;
}

sub get_description {
	my $self = shift;
	return $self->{description};
}

sub get_transcript {
	my $self = shift;
	my ($transcript_name) = @_;
	return $self->{transcript}->{$transcript_name};
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
	my ($start_coord) = @_;
	$self->{start} = $start_coord;
}

sub set_end {
	my $self = shift;
	my ($end_coord) = @_;
	$self->{end} = $end_coord;
}

sub get_chrom {
	my $self = shift;
	return $self->{chrom};
}

sub set_chrom {
	my $self = shift;
	my ($chrom) = @_;
	$self->{chrom} = $chrom;
	# TODO
	#$self->{chrom_index}->{$chrom}->{$gene_id} =
	#					$gene_temp_feature;
}

sub get_strand {
	my $self = shift;
	return $self->{strand};
}

sub set_strand {
	my $self = shift;
	my ($strand) = @_;
	$self->{strand} = $strand;
}

sub get_transcripts_hash {
	my $self = shift;
	$self->{transcript};
}

sub get_transcripts_array {
	my $self = shift;
	return values %{ $self->{transcript} };
}


sub adjust_boundaries_transcript_based {
	my $self = shift;

	$self->set_start( $self->get_transcript_span_start() );
	$self->set_end( $self->get_transcript_span_end() );

}

sub get_transcript_span_start {
	my $self              = shift;
	my $lower_start_value = -1;
	foreach my $ref_transcript ( values %{ $self->{transcript} } ) {
		my $transcript_start = $ref_transcript->get_start();
		$lower_start_value = $transcript_start if ( $transcript_start < $lower_start_value || $lower_start_value == -1 );
	}
	return $lower_start_value;
}

sub get_transcript_span_end {
	my $self             = shift;
	my $higher_end_value = 0;
	foreach my $ref_transcript ( values %{ $self->{transcript} } ) {
		my $transcript_end = $ref_transcript->get_end();
		$higher_end_value = $transcript_end if $transcript_end > $higher_end_value;
	}
	return $higher_end_value;
}

sub toGFF {
	my $self = shift;
	my $str =
	    $self->{chrom} . "\t."
	  . "\tgene\t"
	  . $self->{start} . "\t"
	  . $self->{end} . "\t.\t"
	  . $self->{strand} . "\t.\t" . "ID="
	  . $self->{id} . ";" . "Name="
	  . $self->{name} . ";";

	# Additional attributes
	foreach my $attribName ( keys %{ $self->{attrib} } ) {
		next if( $attribName eq "Name" );
		next if( $attribName eq "ID" );
		$str .= "$attribName=" . $self->{attrib}{$attribName} . ";";
	}
	$str =~ s/$/\n/;

	foreach my $ref_transcript ( values %{ $self->{transcript} } ) {
		$str .= $ref_transcript->toGFF( $self->{chrom}, $self->{id}, $self->{strand} );
	}
	return $str;
}


sub toGTF {
	my $self = shift;
	my $str;
	foreach my $ref_transcript ( values %{ $self->{transcript} } ) {
		$str .= $ref_transcript->toGTF( $self->{chrom}, $self->{id}, $self->{strand} );
	}
	return $str;
}

sub get_num_transcripts {
	my $self = shift;	
	return scalar( values %{ $self->{transcript} } );
}


sub set_attribute {
	my $self = shift;
	my ( $key, $value ) = @_;

	return $self->{attrib}{$key} = $value;
}

sub check_attribute {
	my $self = shift;
	my ($key) = @_;
	return defined $self->{attrib}{$key};
}

sub get_attribute {
	my $self = shift;
	my ($key) = @_;
	if( not defined $self->{attrib}{$key} ){
		die "Gene attribute $key does not exist for gene " 
		. $self->get_id() . " on file " . $self->get_filename() . "\n";
	}
	return $self->{attrib}{$key};
}

sub num_exons {
	my $self = shift;

	my $num_exons = 0;
	foreach my $ref_transcript ( values %{ $self->{transcript} } ) {
		$num_exons += $ref_transcript->num_exons();
	}

	return $num_exons;
}

return 1;
