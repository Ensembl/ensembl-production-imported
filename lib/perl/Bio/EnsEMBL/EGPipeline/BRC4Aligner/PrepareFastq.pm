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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::PrepareFastq;

use strict;
use warnings;
use File::Spec::Functions;

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub run {
  my ($self) = @_;
  
  my $work_dir = $self->param_required('work_dir');
  my $sample_name = $self->param_required("sample_name");

  my $align_meta = $self->param_required('aggregated_aligner_metadata');
  
  my ($file1, $file2);
  if ($align_meta->{is_paired}) {
    $file1 = catfile($work_dir, $sample_name . "_all_1.fastq.gz");
    $file2 = catfile($work_dir, $sample_name . "_all_2.fastq.gz");
  } else {
    $file1 = catfile($work_dir, $sample_name . "_all.fastq.gz");
  }
  
  my $sample_data = {
    sample_seq_file_1 => $file1,
    sample_seq_file_2 => $file2,
  };
  $self->dataflow_output_id($sample_data, 2);
}

1;
