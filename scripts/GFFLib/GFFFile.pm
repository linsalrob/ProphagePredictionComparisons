package GFFFile;
use strict 'vars';
use strict 'refs';

use GFFGene;
use Data::Dumper;
use Carp;

#use Carp::Always;

$Carp::Verbose = 1;

#############
# Constants
my $debug = 0;

my $GENE       = "gene";
my $PSEUDOGENE = "pseudogene";
my $MUORF      = "uORF";

# Transcript types
my $TRANSCRIPT  = "transcript";
my $MRNA        = "mRNA";
my $NCRNA       = "ncRNA";
my $RRNA        = "rRNA";
my $TRNA        = "tRNA";
my $PSEUDO_TRNA = "pseudotRNA";
my $SNRNA       = "snRNA";
my $SCRNA       = "scRNA_encoding";
my $SNORNA      = "snoRNA";
my $MIRNA       = "miRNA";
my $MISCRNA     = "misc_RNA";
my $EXON        = "exon";

my $NON_CODING_EXON     			= "noncoding_exon";
my $CDS  = "CDS";
my $UTR5 = "five_prime_UTR";
my $UTR3 = "three_prime_UTR";

# Genomeview and Ergatis do not like individual IDs for CDS.
# According to those applications those features should share the same ID if
# they belong to the same transcript.
# If the constant below is set to 1, this directive will be imposed when assigng IDs to
# those features.
my $individual_CDS_ids = 1;

sub new {
	my ($filename) = @_;
	my %chroms;
	my $self = { filename => $filename,
				 chroms => %chroms };
	
	bless $self, GFFFile;
	
	return $self;
}

sub get_header_comments{
	my $self = shift;
	return $self->{header_comments};	
}

sub add_exon_to_index {
	my $self = shift;
	
	my ($exon_id, $transcript_id, $gene_id ) = @_;
	
	$self->{exon_index}->{$exon_id} = {
		transcript => $transcript_id,
		gene       => $gene_id
	};
}

sub get_filename {
	my $self = shift;
	return $self->{filename};
}

sub get_chrom_names {
	my $self = shift;
	return keys (%{$self->{chroms}});
}


sub read {
	my $self = shift;
	
	# Add missing features:
	# - add missing genes, transcripts and exons from existing CDS
	# - genes exons or CDS (NOT IMPLEMENTED)
	# - transcripts from exon or CDS (NOT IMPLEMENTED)
	# - exons from CDS (NOT IMPLEMENTED)
	
	my ( $discard_additional_attributes, $add_missing ) = @_;
		
	print STDERR "Reading genes...\n";
	$self->read_genes( $discard_additional_attributes, $add_missing );
	print STDERR "Reading transcripts...\n";
	$self->read_transcripts( $discard_additional_attributes, $add_missing );
	print STDERR "Reading exons...\n";
	$self->read_exons( $discard_additional_attributes, $add_missing );
	print STDERR "Reading UTRs...\n";
	$self->read_utrs( $discard_additional_attributes, $add_missing );
	print STDERR "Reading CDS...\n";
	$self->read_cds( $discard_additional_attributes, $add_missing );
	print STDERR "Post-processing...\n";
	$self->post_processing($add_missing);
}

sub read_genes {
	my $self = shift;
	
	my ( $discard_additional_attributes, $add_missing ) = @_;
	
	open GFFin, $self->{filename}
	or croak "Unable to open file \'" . $self->{filename} . "\' for reading\n";
	
	# Reading only genes
	my $still_on_header = 1;
	while (<GFFin>) {
		my $line = $_;
		
		if ( $line =~ /^\#\#FASTA/ ){
			last;
		} 
		
		if ( $line =~ /^\#/ && $still_on_header){
			$self->{header_comments} .= $line; 
			next;
		}
		next if ( $line =~ /^\#/ || $line =~ /^\s+$/ );
		$still_on_header = 0;
		
		my (
			$chrom,  $feature,     $start, $end,
			$strand, $phase,       $id,    $parent,
			$name,   $description, $attribRef
			) = process_line($line);
			
			
			$self->{chroms}{$chrom} = 1;
			
			my %attr = %{$attribRef};
			
			# Reference to the new feature
			# Used in the iteration to innitialize values and call methods
			# shared by all different type of features
			my $temp_feature;
			
			print STDERR "Line describes feature: $feature\n" if $debug;
			if (   $feature eq $GENE
				|| $feature eq $PSEUDOGENE
				|| $feature eq $MUORF )
			{
				print STDERR "Processing line:\n\t$line\n" if $debug;
				print STDERR "Gene: $id\n"                 if $debug;
				getc()                                     if $debug;
				
				$name = $id if $name eq "";
				
				$temp_feature =
				GFFGene::new( $id, $name, $chrom, $start, $end, $strand, $self );
				$temp_feature->set_filename( $self->{filename} );
				$temp_feature->set_description($description);
				
				$self->{genes}->{$id} = $temp_feature;
				$self->{chrom_index}->{$chrom}->{$id} = $temp_feature;
				
				#######################################
				## General initialization of features
				
				# Set additional attributes
				if ( $discard_additional_attributes == 0 )
				{
					foreach my $currKey ( keys %attr ) {
						$temp_feature->set_attribute( $currKey, $attr{$currKey} );
					}
				}
			}
	}
	close(GFFin);
}

sub read_transcripts {
	
	my $self = shift;
	
	my ( $discard_additional_attributes, $add_missing ) = @_;
	
	# Reading only transcripts
	open GFFin, $self->{filename}
	or croak "Impossible to open file " . $self->{filename} . " for read\n";
	while (<GFFin>) {
		my $line = $_;
		
		if ( $line =~ /^\#\#FASTA/ ){
			last;
		} 
		
		next if ( $line =~ /^\#/ || $line =~ /^\s+$/ );
		
		my (
			$chrom,  $feature,     $start, $end,
			$strand, $phase,       $id,    $parent,
			$name,   $description, $attribRef
			) = process_line($line);
			
			$self->{chroms}{$chrom} = 1;			
			my %attr = %{$attribRef};
			
			# Reference to the new feature
			# Used in the iteration to innitialize values and call methods
			# shared by all different type of features
			my $temp_feature;
			
			if (   $feature eq $TRANSCRIPT
				|| $feature eq $MRNA
				|| $feature eq $RRNA
				|| $feature eq $TRNA
				|| $feature eq $SCRNA
				|| $feature eq $SNRNA
				|| $feature eq $SNORNA
				|| $feature eq $MIRNA
				|| $feature eq $MISCRNA
				|| $feature eq $NCRNA
				|| $feature eq $PSEUDO_TRNA )
			{
				
				print STDERR "Processing line:\n\t$line\n"      if $debug;
				print STDERR "Transcript: $id\tGene: $parent\n" if $debug;
				getc()                                          if $debug;
				
				
				if ( $parent eq '' || not defined $self->{genes}->{$parent} ) {
					if ( $add_missing != 0 ) {
						
						my $gene_name = $id . "-G";
						my $gene_id   = $id . "-G";
						
						my $gene_temp_feature =
						GFFGene::new( $gene_id, $gene_name, $chrom, $start, $end,
							$strand, $self );
						$gene_temp_feature->set_filename( $self->{filename} );
						$self->{genes}->{$gene_id} = $gene_temp_feature;
						$self->{chrom_index}->{$chrom}->{$gene_id} =
						$gene_temp_feature;
						
				
				if ( $discard_additional_attributes == 0 )
				{
					foreach my $currKey ( keys %attr ) {
						$gene_temp_feature->set_attribute( $currKey, $attr{$currKey} );
					}
				}						
						
						$parent = $gene_id;
					}else{
						if( $parent eq '' ){
							croak "Error on file \'"
							. $self->{filename}
							. "\' : The transcript $id does not seem to have a parent\n";
						}else{
							croak "Error on file \'"
							. $self->{filename}
							. "\' : Unable to find parent of $id, looking for $parent\n";
							
							
						}
					}
				}
				
				$temp_feature =
				$self->{genes}->{$parent}
				->add_transcript( $id, $feature, $start, $end );
				$self->{transcript_index}->{$id} = $parent;
				
				#######################################
				## General initialization of features
				# Set additional attributes
				
				#Set additional attributes
				#print ">>>>$discard_additional_attributes<<\n";
				#getc();
							
				if ( $discard_additional_attributes == 0 )
				{
					foreach my $currKey ( keys %attr ) {
						$temp_feature->set_attribute( $currKey, $attr{$currKey} );						
					}
				}
				
			}
			
	}
	close(GFFin);
}

sub read_exons {
	my $self = shift;
	my ( $discard_additional_attributes, $add_missing ) = @_;
	
	# Reading exons
	open GFFin, $self->{filename}
	or croak "Impossible to open file " . $self->{filename} . " for read\n";
	while (<GFFin>) {
		my $line = $_;
		
		if ( $line =~ /^\#\#FASTA/ ){
			last;
		} 		
		
		next if ( $line =~ /^\#/ || $line =~ /^\s+$/ );
		
		print STDERR "Processing line:\n\t$line\n" if $debug;
		my (
			$chrom,  $feature,     $start, $end,
			$strand, $phase,       $id,    $parent,
			$name,   $description, $attribRef
			) = process_line($line);
			
			$self->{chroms}{$chrom} = 1;			
			my %attr = %{$attribRef};
			
			# Reference to the new feature
			# Used in the iteration to innitialize values and call methods
			# shared by all different type of features
			my $temp_feature;
			
			next
			if ( $feature ne $EXON && $feature ne $NON_CODING_EXON);
			
			# The standard GFF format accepts multiple exon Parents, if separated by comma
			my @parents = split ",", $parent;
			
			foreach my $curr_parent (@parents) {
				
				my $curr_transcript = $curr_parent;
				my $curr_gene_name =
				$self->get_gene_name_by_transcript_name($curr_transcript);
				
				if ( not defined($curr_gene_name) ) {
					
					if ( not $add_missing ) {
						croak "Unable to find the transcript and gene associated to"
						. " the exon ID=$id NAME=$name using as parent PARENT=$curr_parent\n";
					}
					else {
						
						my $temp_gene = $self->get_gene($curr_parent);
						
						# If does not have a transcript but has a gene
						if ( defined($temp_gene) ) {
							print STDERR
							"Adding missing transcript to an oprhan exon \'$id\' associate to gene \'$curr_parent\'\n";
							
							my $gene_id = $curr_parent;
							
							# Create new transcript
							my $transcript_id     = $gene_id . "-T";
							my $temp_gene_feature = $self->{genes}->{$gene_id};
							
							$temp_gene_feature->add_transcript(
								$transcript_id, "mRNA",
								$temp_gene_feature->get_start(),
								$temp_gene_feature->get_end()
								);
							$self->{transcript_index}->{$transcript_id} = $gene_id;
							
							$curr_gene_name  = $gene_id;
							$curr_transcript = $transcript_id;
							
						}
						else {
							print STDERR
							"Adding missing gene, transcript to orphan exon \'$id\'\n";
							
							#Generate a random gene id
							my $gene_id = "tmp_gene_" . int( rand(10000) );
							while ( defined $self->{genes}->{$gene_id} ) {
								$gene_id = "tmp_gene_" . int( rand(10000) );
							}
							
							my $temp_gene_feature =
							GFFGene::new( $gene_id, $gene_id, $chrom, $start,
								$end, $strand, $self );
							$temp_gene_feature->set_filename( $self->{filename} );
							$self->{genes}->{$gene_id} = $temp_gene_feature;
							$temp_gene_feature->set_attribute( "Alias",
								$attr{Alias} );
							
							# Transfering attribute from EXON to GENE
							foreach my $currKey ( keys %attr ) {
								$temp_gene_feature->set_attribute( $currKey,
									$attr{$currKey} );
							}
							
							$self->{chrom_index}->{$chrom}->{$id} =
							$temp_gene_feature;
							
							# Create new transcript
							my $transcript_id = $gene_id . "-T";
							$self->{genes}->{$gene_id}
							->add_transcript( $transcript_id,
								"mRNA", $start, $end );
							$self->{transcript_index}->{$transcript_id} = $gene_id;
							
							$curr_gene_name  = $gene_id;
							$curr_transcript = $transcript_id;
							
						}
					}
				}
				
				# If exon doesn't have a name create one or $force_idchange_transcript_exon_cds == 1
				if ( $id eq "" ) {
					$id =
					$parent . "-E"
					. ( $self->{genes}->{$curr_gene_name}
						->get_transcript($curr_transcript)->num_exons() + 1 );
				}
				print STDERR
				"Exon: $id\tTranscript: $curr_transcript\tGene: $curr_gene_name\n"
				if $debug;
				
				$temp_feature =
				$self->{genes}->{$curr_gene_name}
				->get_transcript($curr_transcript)->add_exon( $id, $start, $end );
				$self->{exon_index}->{$id} = {
					transcript => $curr_transcript,
					gene       => $self->{transcript_index}->{$curr_transcript}
				};
			}
			
			#######################################
			## General initialization of features
			
			#	print "Setting attributes >>>$discard_additional_attributes ....\n";
			#	getc();
			
			# Set additional attributes
			#if( $discard_additional_attributes == 0 ){
			foreach my $currKey ( keys %attr ) {
				$temp_feature->set_attribute( $currKey, $attr{$currKey} );
			}
			
			#	print "Setting attributes ....\n";
			#	getc();
			#}
			getc() if $debug;
			
	}
	
	close(GFFin);
	
}

sub read_utrs {
	my $self = shift;
	my ( $discard_additional_attributes, $add_missing ) = @_;
	
	# Reading exons and CDS
	open GFFin, $self->{filename}
	or croak "Impossible to open file " . $self->{filename} . " for read\n";
	while (<GFFin>) {
		my $line = $_;
		
		if ( $line =~ /^\#\#FASTA/ ){
			last;
		} 
		
		
		next if ( $line =~ /^\#/ || $line =~ /^\s+$/ );
		
		print STDERR "Processing line:\n\t$line\n" if $debug;
		my (
			$chrom,  $feature,     $start, $end,
			$strand, $phase,       $id,    $parent,
			$name,   $description, $attribRef
			) = process_line($line);
			
			$self->{chroms}{$chrom} = 1;			
			my %attr = %{$attribRef};
			
			# Reference to the new feature
			# Used in the iteration to innitialize values and call methods
			# shared by all different type of features
			my $temp_feature;
			
			next
			if ( uc($feature) ne uc($UTR5)
				&& uc($feature) ne uc($UTR3) );
			
			if ( uc($feature) eq uc($UTR5) ) {
				my $curr_gene_name =
				$self->get_gene_name_by_transcript_name($parent);
				
				# If it doesn't have a name create one
				if ( $id eq "" ) {
					$id =
					$parent . "-UTR"
					. ( $self->{genes}->{$curr_gene_name}->get_transcript($parent)
						->num_UTRs() + 1 );
				}
				print STDERR
				"UTR5: $id\tTranscript: $parent\tGene: $curr_gene_name\n"
				if $debug;
				
				$temp_feature =
				$self->{genes}->{$curr_gene_name}->get_transcript($parent)
				->add_5UTR( $id, $start, $end );
			}
			
			if ( uc($feature) eq uc($UTR3) ) {
				my $curr_gene_name =
				$self->get_gene_name_by_transcript_name($parent);
				
				# If it doesn't have a name create one
				if ( $id eq "" ) {
					$id =
					$parent . "-UTR"
					. ( $self->{genes}->{$curr_gene_name}->get_transcript($parent)
						->num_UTRs() + 1 );
				}
				print STDERR
				"UTR3: $id\tTranscript: $parent\tGene: $curr_gene_name\n"
				if $debug;
				
				$temp_feature =
				$self->{genes}->{$curr_gene_name}->get_transcript($parent)
				->add_3UTR( $id, $start, $end );
			}
			
			#######################################
			## General initialization of features
			
			#	print "Setting attributes >>>$discard_additional_attributes ....\n";
			#	getc();
			
			# Set additional attributes
			#if( not defined( $discard_additional_attributes ) ||
			#    ( defined( $discard_additional_attributes ) && $discard_additional_attributes == 0 ) ){
			foreach my $currKey ( keys %attr ) {
				$temp_feature->set_attribute( $currKey, $attr{$currKey} );
			}
			
			#	print "Setting attributes ....\n";
			#	getc();
			#}
			getc() if $debug;
	}
	
	close(GFFin);
	
}

sub read_cds {
	my $self = shift;
	my ( $discard_additional_attributes, $add_missing ) = @_;
	
	# Reading exons and CDS
	open GFFin, $self->{filename}
	or croak "Impossible to open file " . $self->{filename} . " for read\n";
	while (<GFFin>) {
		my $line = $_;
		
		if ( $line =~ /^\#\#FASTA/ ){
			last;
		} 		
		
		next if ( $line =~ /^\#/ || $line =~ /^\s+$/ );
		
		print STDERR "Processing line:\n\t$line\n" if $debug;
		my (
			$chrom,  $feature,     $start, $end,
			$strand, $phase,       $id,    $parent,
			$name,   $description, $attribRef
			) = process_line($line);
			
			$self->{chroms}{$chrom} = 1;			
			my %attr = %{$attribRef};
			
			# Reference to the new feature
			# Used in the iteration to innitialize values and call methods
			# shared by all different type of features
			my $temp_feature;
			
			next
			if ( $feature ne $CDS );
			
			my $curr_gene_name = $self->get_gene_name_by_transcript_name($parent);
			
			if ( not defined($curr_gene_name) ) {
				
				if ( not $add_missing ) {
					croak "Unable to find the transcript and gene associated to"
					. " the CDS ID=$id NAME=$name using as parent PARENT=$parent\n";
				}
				else {
					
					my $temp_gene = $self->get_gene($parent);
					
					# If does not have a transcript but has a gene
					if ( defined($temp_gene) ) {
						print STDERR
						"Adding missing transcript to an oprhan CDS \'$id\' associate to gene \'$parent\'\n";
						
						my $gene_id = $parent;
						
						# Check if this gene already have an transcript
						# if it has more than one transcript than DIE, impossible to
						# define which transcript the exon will be associated too
						
						my $temp_gene_feature = $self->{genes}->{$gene_id};
						
						my $num_transcripts = $temp_gene_feature->get_num_transcripts();
						
						if ( $num_transcripts == 1 ) {
							my @trans_arr =
							$temp_gene_feature->get_transcripts_array();
							my $temp_transcript = $trans_arr[0];
							
							croak
							"Undefined value for the first transcript of gene $gene_id\n"
							if not defined($temp_transcript);
							
							$curr_gene_name = $gene_id;
							$parent         = $temp_transcript->get_id();
							
							# If there is no transcripts than create one
						}
						elsif ( $num_transcripts == 0 ) {
							
							# Create new transcript
							my $transcript_id = $gene_id . "-T";
							
							my $temp_transcript =
							$temp_gene_feature->add_transcript(
								$transcript_id, "mRNA",
								$temp_gene_feature->get_start(),
								$temp_gene_feature->get_end()
								);
							$self->{transcript_index}->{$transcript_id} = $gene_id;
							
							$curr_gene_name = $gene_id;
							$parent         = $transcript_id;
						}
						else {
							croak
							"ERROR: Attempting to associate orphan CDS $id to a transcript of gene $gene_id.\n"
							. "But this genes has multiple transcripts (" . $num_transcripts . " transcripts), impossible to define which one will be\n"
							. "associated to the orphan CDS\n";
						}
						
					}
					elsif ( not defined($temp_gene) ) {
						
						# Assuming that CDS represents the whole gene, since there are no gene
						# associated with it.
						# Add missing gene, transcript and EXON
						
						print STDERR
						"Adding missing gene, transcript and exon based on orphan CDS \'$id\'\n";
						
						#Generate a random gene id
						my $gene_id = "tmp_gene_" . int( rand(10000) );
						while ( defined $self->{genes}->{$gene_id} ) {
							$gene_id = "tmp_gene_" . int( rand(10000) );
						}
						
						#Generate a gene id based on CDS id
						#my $gene_id = "$id:$start:$end-gene";
						
						my $temp_gene_feature =
						GFFGene::new( $gene_id, $gene_id, $chrom, $start, $end,
							$strand, $self );
						$temp_gene_feature->set_filename( $self->{filename} );
						$self->{genes}->{$gene_id} = $temp_gene_feature;
						$temp_gene_feature->set_attribute( "Alias", $attr{Alias} );
						
						# Transfering attribute from CDS to GENE
						foreach my $currKey ( keys %attr ) {
							$temp_gene_feature->set_attribute( $currKey,
								$attr{$currKey} );
						}
						
						$self->{chrom_index}->{$chrom}->{$id} = $temp_gene_feature;
						
						# Create new transcript
						my $transcript_id = $gene_id . "-T";
						$self->{genes}->{$gene_id}
						->add_transcript( $transcript_id, "mRNA", $start, $end );
						$self->{transcript_index}->{$transcript_id} = $gene_id;
						
						# Create new exon
						my $exon_id = $gene_id . "-E";
						$self->{genes}->{$gene_id}->get_transcript($transcript_id)
						->add_exon( $exon_id, $start, $end );
						$self->{exon_index}->{$exon_id} = {
							transcript => $transcript_id,
							gene       => $gene_id
						};
						
						$parent         = $transcript_id;
						$curr_gene_name = $gene_id;
					}
				}
			}
			
			# If it doesn't have a name create one
			if ( $id eq "" ) {
				my $serial_num = "";
				$serial_num =
				( $self->{genes}->{$curr_gene_name}->get_transcript($parent)
					->num_CDSs() + 1 )
				if not $individual_CDS_ids;
				$id = $parent . "-P" . $serial_num;
			}
			print STDERR "CDS: $id\tTranscript: $parent\tGene: $curr_gene_name\n"
			if $debug;
			
			$temp_feature =
			$self->{genes}->{$curr_gene_name}->get_transcript($parent)
			->add_CDS( $id, $start, $end, $phase );
			
			
			#######################################
			## General initialization of features
			
			#	print "Setting attributes >>>$discard_additional_attributes ....\n";
			#	getc();
			
			# Set additional attributes
			#if( not defined( $discard_additional_attributes ) ||
			#    ( defined( $discard_additional_attributes ) && $discard_additional_attributes == 0 ) ){
			foreach my $currKey ( keys %attr ) {
				$temp_feature->set_attribute( $currKey, $attr{$currKey} );
			}
			
			#	print "Setting attributes ....\n";
			#	getc();
			#}
			getc() if $debug;
	}
	
	close(GFFin);
	
}

sub toGFF {
	my $self = shift;
	my $str  = "";
	for my $currGene ( values %{ $self->{genes} } ) {
		$str .= $currGene->toGFF();
	}
	return $str;
}

sub write {
	my $self = shift;
	my ($filename) = @_;
	
	open GFFout, ">$filename"
	or croak "Impossible to open file " . $filename . " for write\n";
	
	print GFFout $self->toGFF();
	
	close(GFFout);
}

sub post_processing {
	my $self = shift;
	my ($add_missing) = @_;
	for my $currGene ( values %{ $self->{genes} } ) {
		
		if ( not defined( $currGene->get_transcripts_hash() ) ) {
			
			if ( $add_missing == 0 ) {
				warn "Warning: Gene "
				. $currGene->get_id()
				. " has no transcript. It will be removed!\n";
				delete( $self->{genes}{ $currGene->get_id() } );
				next;
				
			}
			else {
				
				# Create new transcript
				my $gene_id = $currGene->get_id();
				my $start   = $currGene->get_start();
				my $end     = $currGene->get_end();
				
				my $transcript_id = $gene_id . "-T";
				$self->{genes}->{$gene_id}
				->add_transcript( $transcript_id, "mRNA", $start, $end );
				$self->{transcript_index}->{$transcript_id} = $gene_id;
				
				# Create new exon
				my $exon_id = $gene_id . "-E";
				$self->{genes}->{$gene_id}->get_transcript($transcript_id)
				->add_exon( $exon_id, $start, $end );
				$self->{exon_index}->{$exon_id} = {
					transcript => $transcript_id,
					gene       => $gene_id
				};
				
			}
		}
		
		for my $currTrans ( values %{ $currGene->get_transcripts_hash() } ) {
			$currTrans->post_processing($add_missing);
		}
	}
}

sub get_genes_hash {
	my $self = shift;
	return $self->{genes};
}

sub get_genes_hash_by_chrom {
	my $self = shift;
	my ($chrom) = @_;
	return $self->{chrom_index}->{$chrom};
}

sub get_gene {
	my $self = shift;
	my ($gene_name) = @_;
	return $self->{genes}->{$gene_name};
}

sub get_gene_name_by_exon_name {
	my $self = shift;
	my ($exon_name) = @_;
	return $self->{exon_index}->{$exon_name}{gene};
}

sub get_gene_name_by_transcript_name {
	my $self = shift;
	my ($transcript_name) = @_;
	return $self->{transcript_index}->{$transcript_name};
}

sub get_gene_by_exon_name {
	my $self = shift;
	my ($exon_name) = @_;
	
	my $gene_name = $self->{exon_index}->{$exon_name}{gene};
	return $self->{genes}->{$gene_name};
}

sub get_gene_by_transcript_name {
	my $self = shift;
	my ($transcript_name) = @_;
	
	my $gene_name = $self->{transcript_index}->{$transcript_name};
	return $self->{genes}->{$gene_name};
}

sub get_transcript_by_name {
	my $self = shift;
	my ($transcript_name) = @_;
	
	#print ">>>>" . $transcript_name . "<<<<<\n";
	#print Dumper( $self->get_gene_by_transcript_name( $transcript_name ) );
	#getc();
	my $gene = $self->get_gene_by_transcript_name($transcript_name);
	return 0 if ( not defined $gene );
	
	return $gene->get_transcript($transcript_name);
}

sub get_transcript_by_exon_name {
	my $self = shift;
	my ($exon_name) = @_;
	
	#print ">>>>" . $exon_name . "<<<<<\n";
	#print Dumper( $self->{exon_index}->{$exon_name} );
	#getc();
	
	my $gene_name = $self->{exon_index}->{$exon_name}->{gene};
	return undef if not defined($gene_name);
	
	my $transcript_name = $self->{exon_index}->{$exon_name}->{transcript};
	return undef if not defined($transcript_name);
	
	return $self->{genes}->{$gene_name}->get_transcript($transcript_name);
}

return 1;

sub process_line {
	my ($line) = @_;
	
	if (
		(
			my (
				$chrom, $void1,  $feature, $start, $end,
				$void2, $strand, $phase, $info
				)
			= split "\t", $line
			) == 9
		)
	{
		
		my ($id)     = ( $info =~ /ID=([\w\W]+?)[;\n]/ );
		my ($parent) = ( $info =~ /Parent=([\w\W]+?)[;\n]/ );
		my ($name)   = ( $info =~ /Name=([\w\W]+?)[;\n]/ );
		
		chomp($info);
		
		#Retrieving all attributes associated to the feature
		my %attr;
		my @separate_semicolon = split ";", $info;
		foreach my $currAttr (@separate_semicolon) {
			my ( $key, $value );
			if ( ( $key, $value ) = ( $currAttr =~ /\s*([\w\W]+?)=([\w\W]+)/ ) )
			{
				#				if (   $key ne "ID"
				#					&& $key ne "Parent"
				#					&& $key ne "Name" )
				#				{
				$attr{$key} = $value;
				#				}
			}
		}
		
		# Same description stores the gene product
		my ($description) = ( $info =~ /description=([\w\W]+?)[;\n]/ );
		
		return (
			$chrom,  $feature,     $start, $end,
			$strand, $phase,       $id,    $parent,
			$name,   $description, \%attr
			);
	}
	else {
		print STDERR
		"\nThe following line does not have 9 columns separated by <tab>, a GFF format requirement.\n";
		print STDERR
		"Check if the columns are seperated by multiple <spaces> instead.\n";
		croak "$line";
	}
}
