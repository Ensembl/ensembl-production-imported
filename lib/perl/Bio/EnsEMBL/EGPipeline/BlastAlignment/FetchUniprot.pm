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

Bio::EnsEMBL::EGPipeline::BlastAlignment::FetchUniprot

=head1 DESCRIPTION

Download UniProtKB from ftp site and convert into Fasta format.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::BlastAlignment::FetchUniprot;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::FetchExternal');

use Bio::SeqIO;
use File::Path qw(make_path);
use File::Spec::Functions qw(catdir);
use IO::Uncompress::Gunzip qw(gunzip);

sub param_defaults {
  my ($self) = @_;
  return {
    'ebi_path'     => $self->private_conf('ENSEMBL_UNIPROT_EBI_PATH'),
    'ftp_uri'      => 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase',
    'uniprot_db'   => 'sprot',
    'file_varname' => 'uniprot_fasta_file',
  };
}

sub fetch_input {
  my ($self) = @_;
  my $out_dir         = $self->param_required('out_dir');
  my $uniprot_db      = $self->param_required('uniprot_db');
  my $taxonomic_level = $self->param('taxonomic_level');
  
  if (!-e $out_dir) {
    $self->warning("Output directory '$out_dir' does not exist. I shall create it.");
    make_path($out_dir) or $self->throw("Failed to create output directory '$out_dir'");
  }
  
  my ($uniprot_file, $output_file);
  if ($taxonomic_level) {
    $uniprot_file = "uniprot_$uniprot_db\_$taxonomic_level.dat.gz";
    $output_file  = catdir($out_dir, "$uniprot_db\_$taxonomic_level.fa");
  } else {
    $uniprot_file = "uniprot_$uniprot_db.fasta.gz";
    $output_file  = catdir($out_dir, "uniprot_$uniprot_db.fasta");
  }
  
  $self->param('uniprot_file', $uniprot_file);
  $self->param('output_file',  $output_file);
}

sub run {
  my ($self) = @_;
  my $ebi_path        = $self->param_required('ebi_path');
  my $ftp_uri         = $self->param_required('ftp_uri');
  my $out_dir         = $self->param_required('out_dir');
  my $uniprot_file    = $self->param_required('uniprot_file');
  my $output_file     = $self->param_required('output_file');
  my $taxonomic_level = $self->param('taxonomic_level');
  
  if ($taxonomic_level) {
    $ebi_path = catdir($ebi_path, 'taxonomic_divisions');
    $ftp_uri .= '/taxonomic_divisions';
  } else {
    $ebi_path = catdir($ebi_path, 'complete');
    $ftp_uri .= '/complete';
  }
  
  my $ebi_file   = catdir($ebi_path, $uniprot_file);
  my $local_file = catdir($out_dir,  $uniprot_file);
  
  if (-e $ebi_file) {
    $self->fetch_ebi_file($ebi_file, $local_file);
  } else {
    my $ftp = $self->get_ftp($ftp_uri);
    $self->fetch_ftp_file($ftp, $uniprot_file, $local_file);
  }
  
  if ($taxonomic_level) {
    $self->convert_to_fasta($local_file, 'swiss', $output_file);
  } else {
    $self->convert_to_fasta($local_file, 'fasta', $output_file);
  }
}

sub convert_to_fasta {
  my ($self, $file, $format, $output_file) = @_;
  
  my $seq_out = Bio::SeqIO->new(
    -file   => '>'.$output_file,
    -format => 'fasta',
  );
  
  # Don't forget these are gzipped files...
  open(F, "gunzip -c $file |");
  my $seq_in = Bio::SeqIO->new(
    -fh     => \*F,
    -format => $format,
  );
  
  while (my $inseq = $seq_in->next_seq) {
    my (undef, $primary_id, undef) = split(/\|/, $inseq->display_id);
    $inseq->display_id($primary_id);
    
    $seq_out->write_seq($inseq);
  }
}

1;
