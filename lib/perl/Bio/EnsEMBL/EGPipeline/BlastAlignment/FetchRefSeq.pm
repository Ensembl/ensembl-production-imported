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

Bio::EnsEMBL::EGPipeline::BlastAlignment::FetchRefSeq

=head1 DESCRIPTION

Download RefSeq files from ftp site and consolidate into a single fasta file.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::BlastAlignment::FetchRefSeq;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::FetchExternal');

use File::Basename qw(fileparse);
use File::Path qw(make_path);
use File::Spec::Functions qw(catdir);
use IO::Uncompress::Gunzip qw(gunzip);
use Path::Tiny qw(path);

sub param_defaults {
  my ($self) = @_;
  return {
    'ebi_path'     => $self->private_conf('ENSEMBL_REFSEQ_EBI_PATH'),
    'ftp_uri'      => 'ftp://ftp.ncbi.nlm.nih.gov/refseq/release',
    'data_type'    => 'pep',
    'file_varname' => 'refseq_fasta_file',
  };
}

sub fetch_input {
  my ($self) = @_;
  my $out_dir         = $self->param_required('out_dir');
  my $taxonomic_level = $self->param_required('taxonomic_level');
  my $data_type       = $self->param_required('data_type');
  
  if (!-e $out_dir) {
    $self->warning("Output directory '$out_dir' does not exist. I shall create it.");
    make_path($out_dir) or $self->throw("Failed to create output directory '$out_dir'");
  }
  
  $self->param('ebi_path', catdir($self->param('ebi_path'), $taxonomic_level));
  $self->param('ftp_uri',  $self->param('ftp_uri') . "/$taxonomic_level");
  
  my $refseq_files;
  if ($data_type eq 'pep') {
    $refseq_files = "$taxonomic_level.*.protein.faa.gz";
  } else {
    $refseq_files = "$taxonomic_level.*.rna.fna.gz";
  }
  my $output_file = catdir($out_dir, "refseq_$taxonomic_level\_$data_type.fa");
  
  $self->param('refseq_files', $refseq_files);
  $self->param('output_file', $output_file);
}

sub run {
  my ($self) = @_;
  my $ebi_path     = $self->param_required('ebi_path');
  my $ftp_uri      = $self->param_required('ftp_uri');
  my $out_dir      = $self->param_required('out_dir');
  my $refseq_files = $self->param_required('refseq_files');
  my $output_file  = $self->param_required('output_file');
  
  my $ebi_files = catdir($ebi_path, $refseq_files);
  
  # The trickiness here is that the refseq data is stored in multiple
  # files. We will iterate over them, appending to the output file;
  # so need to make sure any existing file is removed.
  if (-e $output_file) {
    unlink $output_file;
  }
  
  my @ebi_files = glob($ebi_files);
  my @local_files;
  
  if (scalar(@ebi_files)) {
    foreach my $ebi_file (@ebi_files) {
      my $filename = fileparse($ebi_file);
      my $local_file = catdir($out_dir, $filename);
      
      $self->fetch_ebi_file($ebi_file, $local_file);
      push @local_files, $local_file;
    }
  } else {
    my $ftp = $self->get_ftp($ftp_uri);
    @local_files = @{ $self->fetch_ftp_files($ftp, $refseq_files, $out_dir) };
  }
  
  my $output = path($output_file);
  foreach my $local_file (@local_files) {
    gunzip $local_file => "$local_file.tmp";
    
    my $data = path("$local_file.tmp")->slurp;
    $output->append($data);
    
    unlink "$local_file.tmp";
  }
}

1;
