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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::AlignSequence;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::BRC4Aligner::Base');

sub param_defaults {
  my ($self) = @_;
  
  return {
    'threads'  => 4,
    'run_mode' => 'default',
    'gtf_file' => undef,
    'clean_up' => 1,
    'store_cmd' => 1,
    dnaseq      => 0,
  };
}

sub fetch_input {
  my ($self) = @_;
  
  if (defined $self->param('escape_branch') and $self->input_job->retry_count > 0 and $self->input_job->retry_count >= $self->input_job->analysis->max_retry_count) {
    $self->dataflow_output_id($self->input_id, $self->param('escape_branch'));
    $self->input_job->autoflow(0);
    $self->complete_early("Failure probably due to memory limit, retrying with a higher limit.");
  }
  
  my $aligner_class = $self->param_required('aligner_class');
  my $threads       = $self->param_required('threads');
  my $run_mode      = $self->param_required('run_mode');
  my $aligner_dir   = $self->param('aligner_dir');
  my $samtools_dir  = $self->param('samtools_dir');
  my $max_intron    = $self->param('max_intron');
  my $gtf_file      = $self->param('gtf_file');
  my $strandness    = $self->param('strandness');
  my $dnaseq        = $self->param('dnaseq');

  # DNA-Seq special changes
  $strandness = undef if $dnaseq;
  my $no_spliced = $dnaseq ? 1 : 0;
  my $number_primary = $dnaseq ? 1 : undef;
  
  my $max_intron_length = $self->max_intron_length() if $max_intron;
  
  eval "require $aligner_class";
  
  my $aligner_object = $aligner_class->new(
    -aligner_dir       => $aligner_dir,
    -samtools_dir      => $samtools_dir,
    -threads           => $threads,
    -run_mode          => $run_mode,
    -gtf_file          => $gtf_file,
    -max_intron_length => $max_intron_length,
    -strandness        => $strandness,
    -no_spliced        => $no_spliced,
    -number_primary    => $number_primary,
  );
  
  $self->param('aligner_object', $aligner_object);
}

sub run {
  my ($self) = @_;
  
  my $aligner     = $self->param_required('aligner_object');
  my $genome_file = $self->param_required('genome_file');
  my $seq_file_1  = $self->param_required('seq_file_1');
  my $seq_file_2  = $self->param('seq_file_2') || undef;
  my $clean_up    = $self->param_required('clean_up');
  my $sam_file    = $self->param('sam_file');
  
  $sam_file //= "$seq_file_1.sam";
  
  # Align to create a SAM file
  $aligner->align($genome_file, $sam_file, $seq_file_1, $seq_file_2);
  $aligner->dummy(0);
  $self->param('sam_file', $sam_file);
  
  my $index_cmds = $self->param('index_cmds') || [];
  my $align_cmds = $aligner->align_cmds;
  my $version    = $aligner->version;
  
  $self->param('cmds', join("; ", (@$index_cmds, @$align_cmds)));
  $self->param('version', $version);
  $self->param('sam_file', $sam_file);
}

sub write_output {
  my ($self) = @_;
  
  if ($self->param("store_cmd")) {
    my $align_cmds = {
      cmds => $self->param('cmds'),
      sample_name => $self->param('sample_name'),
      'version'  => $self->param('version'),
    };
    $self->store_align_cmds($align_cmds);
  }
  
  my $dataflow_output = {
    'sam_file' => $self->param('sam_file'),
  };
  
  $self->dataflow_output_id($dataflow_output,  2);
}

sub max_intron_length {
  my ($self) = @_;
  
  my $max_intron_length = 0;
  
  my $dba = $self->get_DBAdaptor('core');
  my $ta = $dba->get_adaptor('Transcript');
  
  my $transcripts = $ta->fetch_all_by_biotype('protein_coding');

  my $ceiling = 500_000;
  
  foreach my $transcript (@$transcripts) {
    my $introns = $transcript->get_all_Introns();
    foreach my $intron (@$introns) {
      if ($intron->length() > $max_intron_length) {
        $max_intron_length = $intron->length();
	    }
    }
    if ($max_intron_length > $ceiling) {
      $max_intron_length = $ceiling;
      last;
    }
  }
  $max_intron_length = int(1.5 * $max_intron_length);
  $max_intron_length = $ceiling if $max_intron_length > $ceiling;
  
  return $max_intron_length;
}

1;
