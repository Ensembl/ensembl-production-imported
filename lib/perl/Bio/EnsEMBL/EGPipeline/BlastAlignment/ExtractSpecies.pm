=head1 LICENSE

See the NOTICE file distributed with this work for additional information
regarding copyright ownership.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=pod

=head1 NAME

Bio::EnsEMBL::EGPipeline::BlastAlignment::ExtractSpecies

=head1 DESCRIPTION

Extract entries for a particular species, from Fasta format file.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::BlastAlignment::ExtractSpecies;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::SeqIO;

sub param_defaults {
  return {
    'file_varname' => 'species_fasta_file',
  };
}

sub fetch_input {
  my ($self) = @_;
  my $fasta_file     = $self->param_required('fasta_file');
  my $species        = $self->param('species');
  my $source_species = $self->param('source_species');
 
  if (! defined $source_species) {
    if (! defined $species) {
      $self->throw('-species or -source_species parameter is required');
    } else {
      $source_species = $species;
      $self->param('source_species', $source_species);
    }
  }

  (my $output_file = $fasta_file) =~ s/(\.\w+)$/_$source_species$1/;
  

  my $trinomial_rex = qr/[A-Z]{1}[a-z]+\s[a-z]+\s[a-z0-9]+/;
  my $binomial_rex  = qr/[A-Z]{1}[a-z]+\s[a-z]+/;

  my $trinomial_name;
  my $binomial_name;
  if ($source_species=~m/^([a-z]+_[a-z]+)_[a-z0-9]+$/) {
    $trinomial_name = $source_species;
    $binomial_name  = $1;
    $trinomial_name=~s/_/ /g;
    $trinomial_name = ucfirst($trinomial_name);
  } elsif ($source_species=~m/^([a-z]+_[a-z]+)$/) {
    $binomial_name  = $1;
  } else {
    $self->throw("Naming pattern was not recognised, you must specify a new pattern for the code to work. Source species: '$source_species'");
  }

  $binomial_name=~s/_/ /g;
  $binomial_name  = ucfirst($binomial_name);
  
  # Get taxon id
  my $taxon_id = $self->get_taxon_id($species);
  $self->param('taxon_id', $taxon_id);
  
  $self->param('trinomial_rex', $trinomial_rex);
  $self->param('binomial_rex', $binomial_rex);
  $self->param('trinomial_name', $trinomial_name);
  $self->param('binomial_name', $binomial_name);

  $self->param('output_file', $output_file);
}

sub get_taxon_id {
  my ($self, $species) = @_;
  
  my $dba = $self->core_dba();
  my $sql = 'SELECT meta_value FROM meta WHERE meta_key = ?';
  my $dbh = $dba->dbc->db_handle;
  my $sth = $dbh->prepare($sql);
  $sth->execute("species.taxonomy_id");
  
  my ($taxon_id) = $sth->fetchrow_array;
  
  return $taxon_id;
}

sub run {
  my ($self) = @_;
  my $fasta_file        = $self->param_required('fasta_file');
  my $output_file       = $self->param_required('output_file');
  my $data_source       = $self->param_required('data_source');
  my $data_type         = $self->param_required('data_type');
  
  my $seq_out = Bio::SeqIO->new(
    -file   => '>'.$output_file,
    -format => 'fasta',
  );

  open(F, $fasta_file);
  my $seq_in = Bio::SeqIO->new(
    -fh     => \*F,
    -format => 'fasta',
  );

  while (my $inseq = $seq_in->next_seq) {
    if ($data_source eq 'uniprot') {
      $self->parse_uniprot($seq_out, $inseq);
    } elsif ($data_source eq 'refseq') {
      $self->parse_refseq($data_type, $seq_out, $inseq);
    }
  }
}

sub write_output {
  my ($self) = @_;
  my $output_file  = $self->param('output_file');
  my $file_varname = $self->param_required('file_varname');
  
  if (-s $output_file) {
    $self->dataflow_output_id({$file_varname => $output_file}, 1);
  } else {
    $self->warning("No data for $output_file.");
    $self->input_job->autoflow(0);
  }
}

sub parse_uniprot {
  my ($self, $seq_out, $inseq) = @_;
  my $trinomial_rex  = $self->param('trinomial_rex');
  my $binomial_rex   = $self->param('binomial_rex');
  my $trinomial_name = $self->param('trinomial_name');
  my $binomial_name  = $self->param('binomial_name');
  my $taxon_id       = $self->param('taxon_id');
  
  #warn "$trinomial_name:$binomial_name";
  my $full_desc = $inseq->desc;
  
  # Some 'descriptions' have "OS=<species>" embedded within them,
  # presumably because UniProt records were used as labels, without
  # appropriate filtering.
  my @OSs = $full_desc =~ /OS=(\w+)/g;
  if (scalar(@OSs) == 2) { 
    $full_desc =~ s/\s+OS=.*?SV=\d+//;
  } elsif (scalar(@OSs) > 2) { 
    $self->throw("More than two 'OS=' sections in description: $full_desc");
  }
  
  if($full_desc=~m/OS=$trinomial_rex\sOX/ and $trinomial_name){
  	  if($full_desc=~m/OS=$trinomial_name\sOX/){
  	  	  my ($desc, $version) = $full_desc =~ /^(.*)\s+OS=.*SV=(\d+)/;
  	  	  $inseq->desc(join('|', map { $_ || '' } ($inseq->display_id(), $desc, $version)));
  	  	  $seq_out->write_seq($inseq);
	  }
  }elsif($full_desc=~/OS=$trinomial_rex\sOX/){
  	  if($full_desc=~m/OS=$binomial_name/){
  	  	  #IF you wish descriptions with trinomial names in the heading, 
  	  	  #to be added to species with binomial name, please add code below.
  	  	  
	  }
  }elsif($full_desc=~m/OS=$binomial_rex\sOX/){
  	  if($full_desc=~m/OS=$binomial_name\sOX/){
  	  	  my ($desc, $version) = $full_desc =~ /^(.*)\s+OS=.*SV=(\d+)/;
  	  	  $inseq->desc(join('|', map { $_ || '' } ($inseq->display_id(), $desc, $version)));
  	  	  $seq_out->write_seq($inseq);
	  }
  }elsif($taxon_id and $full_desc=~m/OX=$taxon_id\s/){
    # Use the taxon id, more stable 
    my ($desc, $version) = $full_desc =~ /^(.*)\s+OS=.*SV=(\d+)/;
    $inseq->desc(join('|', map { $_ || '' } ($inseq->display_id(), $desc, $version)));
    $seq_out->write_seq($inseq);
  }
}

sub parse_refseq {
  my ($self, $data_type, $seq_out, $inseq) = @_;
  
  if ($self->species_match($inseq->desc, $data_type)) {
    my ($primary_id, $version) = $inseq->display_id =~ /(\w+)\.(\d+)|[^\|]+$/;
    $inseq->display_id("$primary_id.$version");
    
    my ($desc) = $inseq->desc =~ /^(.*)\s+(?:\[|\().*/;
    $desc =~ s/PREDICTED:\s+//;
    $inseq->desc(join('|', map { $_ || '' } ("$primary_id.$version", $desc, $version)));

    $seq_out->write_seq($inseq);
  }
}

sub species_match {
  my ($self, $desc, $data_type) = @_;

  my $trinomial_rex  = $self->param('trinomial_rex');
  my $binomial_rex   = $self->param('binomial_rex');
  my $trinomial_name = $self->param('trinomial_name');
  my $binomial_name  = $self->param('binomial_name');
  
  if ($data_type eq 'pep') {
	  if($desc=~m/\[$trinomial_rex\]/ and $trinomial_name){
		  if($desc=~m/\[$trinomial_name\]/){
			return 1; 
		  }
	  }elsif($desc=~/\[$trinomial_rex\]/){
		  if($desc=~m/\[$binomial_name\]/){
			return 0; 
		  }
	  }elsif($desc=~m/\[$binomial_rex\]/){
		  if($desc=~m/\[$binomial_name\]/){
			 return 1; 
		  }
	  }
  } else {
    	  if($desc=~m/^$trinomial_rex/ and $trinomial_name){
		  if($desc=~m/^$trinomial_name/){
			return 1; 
		  }
	  }elsif($desc=~/^$trinomial_rex/){
		  if($desc=~m/^$binomial_name/){
			return 0; 
		  }
	  }elsif($desc=~m/^$binomial_rex/){
		  if($desc=~m/^$binomial_name/){
			 return 1; 
		  }
	  }
  }
}

1;
