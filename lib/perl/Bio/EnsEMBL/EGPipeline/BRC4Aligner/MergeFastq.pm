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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::MergeFastq;

use strict;
use warnings;
use File::Spec::Functions;

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub run {
  my ($self) = @_;
  
  my $work_dir = $self->param_required('work_dir');
  my $sample_name = $self->param_required("sample_name");
  my $seq_files_1 = $self->param("run_seq_files_1");
  my $seq_files_2 = $self->param("run_seq_files_2");
  
  die "No seq_file 1 defined" if not defined $seq_files_1;
  
  my ($file1, $file2) = $self->merge_fastq($seq_files_1, $seq_files_2, $work_dir, $sample_name);
  
  my $merged_data = {
    sample_seq_file_1 => $file1,
    sample_seq_file_2 => $file2,
  };
  $self->dataflow_output_id($merged_data, 2);
}

sub merge_fastq {
  my ($self, $files1, $files2, $dir, $sample_name) = @_;
  
  # Remove empty files 2
  @$files2 = grep { $_ } @$files2;
  
  # Checks
  if (not $files1 or scalar(@$files1) == 0) {
    die "No files given to merge!";
  }
  if (scalar(@$files2) > 0 and scalar(@$files1) != scalar(@$files2)) {
    die "Paired files: not the same number of files provided! " . scalar(@$files1) . " vs " . scalar(@$files2) . "(@$files1, @$files2)";
  }
  
  # Prepare the files names
  my ($output1, $output2);
  if (scalar @$files2) {
    $output1 = catfile($dir, $sample_name . "_all_1.fastq.gz");
    $output2 = catfile($dir, $sample_name . "_all_2.fastq.gz");
  } else {
    $output1 = catfile($dir, $sample_name . "_all.fastq.gz");
  }

  # Only one run: no need to merge
  if (scalar(@$files1) == 1) {
    rename $files1->[0], $output1;
    rename $files2->[0], $output2 if @$files2;
  } else {
    
    # Merge the files in order
    $self->cat_files($files1, $output1) if not -s $output1;
    $self->cat_files($files2, $output2) if @$files2 and not -s $output2;
    
    # Delete the merged files
    #unlink @$files1, @$files2;
  }
  
  return $output1, $output2;
}

sub cat_files {
  my ($self, $files, $output) = @_;
  
  return if not $files or @$files == 0;
  return if -s $output;
  
  my $files_str = join(" ", sort @$files);
  
  system("cat $files_str > $output");
}

1;
