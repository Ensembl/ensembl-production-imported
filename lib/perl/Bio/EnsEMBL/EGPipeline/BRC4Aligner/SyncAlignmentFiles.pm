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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::SyncAlignmentFiles;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::EnsEMBL::ENA::SRA::BaseSraAdaptor qw(get_adaptor);

use File::Path qw(make_path);
use File::Spec::Functions qw(catdir catfile);
use File::Copy;
use JSON;

sub fetch_input {
  my ($self) = @_;
  
  my $component = $self->param("component");
  my $organism = $self->param("organism");
  $component = 'component' if not $component;
  $organism = $self->param("species") if not $organism;
  my $dataset = $self->param_required("study_name");
  my $sample = $self->param_required("sample_name");
  
  my $alignment_root_dir = $self->param_required("alignments_dir");
  my $sample_dir = $self->param_required("sample_dir");
  
  # Check that we can get the old alignment files
  my $old_dir = catdir($alignment_root_dir, $component, $organism, $dataset, $sample);
  my $old_metadata_file = catfile($old_dir, 'metadata.json');
  my $old_bam_file = catfile($old_dir, 'results_name.bam');
  
  if (not -e $old_dir) {
    $self->throw("Old alignment dir doesn't exist: $old_dir");
  }
  if (not -s $old_metadata_file) {
    $self->throw("Old alignment metadata doesn't exist: $old_metadata_file");
  }
  if (not -s $old_bam_file) {
    $self->throw("Old alignment bam file doesn't exist: $old_bam_file");
  }
  
  $self->param('old_bam_file', $old_bam_file);
  $self->param('old_metadata_file', $old_metadata_file);
}

sub run {
  my ($self) = @_;
  
  my $sample_dir = $self->param('sample_dir');
  my $old_bam_file = $self->param('old_bam_file');
  my $new_bam_file = $self->param('sorted_bam_file');
  my $old_metadata_file = $self->param('old_metadata_file');
  my $new_metadata_file = $self->param('metadata_file');
  
  if (not -s $sample_dir) {
    make_path($sample_dir);
  }
  if (not -s $new_metadata_file) {
    copy $old_metadata_file, $new_metadata_file;
  }
  if (not -s $new_bam_file) {
    copy $old_bam_file, $new_bam_file;
  }
  
  # Checksum? TODO
}

sub write_output {
  my ($self) = @_;
  $self->dataflow_output_id({}, 2);
}

1;

