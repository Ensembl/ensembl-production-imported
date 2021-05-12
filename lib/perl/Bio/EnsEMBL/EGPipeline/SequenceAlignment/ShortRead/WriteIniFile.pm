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

package Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::WriteIniFile;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use File::Spec::Functions qw(catdir);

sub param_defaults {
  my ($self) = @_;
  
  return {
    'bigwig'   => 0,
    'ini_type' => 'rnaseq_align',
  };
}

sub run {
  my ($self) = @_;
  
  my $species     = $self->param_required('species');
  my $results_dir = $self->param_required('results_dir');
  my $bigwig      = $self->param_required('bigwig');
  my $ini_type    = $self->param_required('ini_type');
  my $merge_level = $self->param_required('merge_level');
  my $merge_id    = $self->param_required('merge_id');
  my $merge_label = $self->param_required('merge_label');
  my $run_ids     = $self->param_required('run_ids');
  
  # Get assembly name
  my $dba = $self->core_dba;
  my $assembly = $dba->get_MetaContainer()->single_value_by_key('assembly.default');
  
  my ($name, $description) = ($merge_label, '');
  my $caption = $ini_type eq 'rnaseq_align' ? 'RNA-seq alignment' : $name;
  
  if (length($run_ids) > 0) {
    ($description) = $self->sra_desc($merge_id, $merge_level, $run_ids);
  }
  
  my ($header, $data_file, $source_type);
  if ($bigwig) {
    $header = '[ENSEMBL_INTERNAL_BIGWIG_SOURCES]';
    $data_file = $self->param('bw_file');
    $source_type = 'bigWig';
  } else {
    $header = '[ENSEMBL_INTERNAL_BAM_SOURCES]';
    $data_file = $self->param('merged_bam_file');
    $source_type = 'bam';
  }
  
  $header .= "\n$merge_id = $ini_type\n";
  my $body = 
    "[$merge_id]\n".
    "source_name = $name\n".
    "caption     = $caption\n".
    "description = $description\n".
    "source_url  = $data_file\n".
    "source_type = $source_type\n".
    "display     = off\n";
  
  my $file_name = sprintf("/%s_%s.ini", $merge_id, $assembly);
  my $ini_file = catdir($results_dir, $file_name);
  open (my $fh, '>', $ini_file) or die "Failed to open file '$ini_file'";
	print $fh "$header\n$body\n";
  close($fh);
  
  $self->param('ini_file', $ini_file);
}

sub write_output {
  my ($self) = @_;
  
  $self->dataflow_output_id({'ini_file' => $self->param('ini_file')}, 1);
}

sub sra_desc {
  my ($self, $merge_id, $merge_level, $run_ids) = @_;
  
  my $description = ucfirst(lc($merge_level)) . ' ' . $self->ena_link($merge_id);
  
  if ($merge_level !~ /^run$/i) {
    my @run_ids = split(/,/, $run_ids);
    my @run_links = map { $self->ena_link($_) } @run_ids;
    
    if (scalar @run_ids == 1) {
      $description .= ' (Run: ';
    } else {
      $description .= ' (Runs: ';
    }
    $description .= join(', ', @run_links).')';
  }
  
  return ($description);
}

sub ena_link {
  my ($self, $id) = @_;
  return "<a href='https://www.ebi.ac.uk/ena/data/view/$id' target='_blank'>$id</a>"
}

1;
