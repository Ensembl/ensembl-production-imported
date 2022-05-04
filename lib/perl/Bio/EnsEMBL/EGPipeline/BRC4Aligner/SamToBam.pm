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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::SamToBam;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::BRC4Aligner::Base');
use Bio::EnsEMBL::EGPipeline::Common::Aligner;

sub param_defaults {
  my ($self) = @_;
  
  return {
    'threads'   => 4,
    'memory'    => 8000,
    'skip_cleanup'  => 0,
  };
}

sub fetch_input {
  my ($self) = @_;
  
  my $samtools_dir = $self->param('samtools_dir');
  my $threads      = $self->param_required('threads');
  my $memory       = $self->param_required('memory');
  
  my $aligner_object = Bio::EnsEMBL::EGPipeline::Common::Aligner->new(
    -samtools_dir => $samtools_dir,
    -threads      => $threads,
    -memory       => $memory,
  );

  $self->param('aligner_object', $aligner_object);
}

sub run {
  my ($self) = @_;
  
  my $aligner  = $self->param_required('aligner_object');
  my $sam_file = $self->param_required('sam_file');
  my $skip_cleanup = $self->param_required('skip_cleanup')
  my $final_bam_file = $self->param('final_bam_file');
  
  # Can we reuse some files?
  if (-s $final_bam_file) {
    warn("Bam file already exists: $final_bam_file");
    $aligner->dummy(1);
  }
  warn("Create bam file: $final_bam_file");
  
  # Convert
  my $bam_file = $aligner->sam_to_bam($sam_file, $final_bam_file);
  unlink $sam_file unless $skip_cleanup;
  $aligner->dummy(0);
  
  my $align_cmds = $aligner->align_cmds;
  
  $self->param('output_bam_file', $bam_file);
  $self->param('cmds', join("; ", @$align_cmds));
}

sub write_output {
  my ($self) = @_;
  
  my $align_cmds = {
    cmds => $self->param('cmds'),
    sample_name => $self->param('sample_name'),
  };
  $self->store_align_cmds($align_cmds);
  
  my $dataflow_output_to_next = {
    'output_bam_file' => $self->param('output_bam_file'),
  };
  
  $self->dataflow_output_id($dataflow_output_to_next,  1);
}

1;
