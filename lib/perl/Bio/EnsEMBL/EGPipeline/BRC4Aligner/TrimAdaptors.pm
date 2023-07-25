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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::TrimAdaptors;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Capture::Tiny ':all';
use Path::Tiny qw(path);
use File::Basename qw(dirname);
use File::Spec::Functions qw(catdir);
use File::Copy qw(move);

sub param_defaults {
  my ($self) = @_;
  
  return {
    'threads' => 1,
  };
}

sub run {
  my ($self) = @_;
  
  my $seq1 = $self->param_required('seq_file_1');
  my $seq2 = $self->param('seq_file_2');

  my $sub_seq1 = $seq1; $sub_seq1 =~ s/(\..+$)/_trim$1/;
  my $sub_seq2 = $seq2; $sub_seq2 =~ s/(\..+$)/_trim$1/;

  my $old_seq1 = $seq1; $old_seq1 =~ s/(\..+$)/_old$1/;
  my $old_seq2 = $seq2; $old_seq2 =~ s/(\..+$)/_old$1/;

  # Trim then replace the file with the trimmed version,
  # and keep the old file
  if ($seq2) {
    $self->trim_reads_paired($seq1, $seq2, $sub_seq1, $sub_seq2);
    move $seq1, $old_seq1;
    move $sub_seq1, $seq1;
    move $seq2, $old_seq2;
    move $sub_seq2, $seq2;
  } else {
    $self->trim_reads_single($seq1, $sub_seq1);
    move $seq1, $old_seq1;
    move $sub_seq1, $seq1;
  }
}

sub write_output {
  my ($self) = @_;
  
  my $output = {
    sample_seq_file_1 => $self->param('seq_file_1'),
    sample_seq_file_2 => $self->param('seq_file_2'),
  };
  $self->dataflow_output_id($output,  2);
}

sub trim_reads_single {
  my ($self, $seq, $seq_trim) = @_;
  return if -s $seq_trim;

  my $trim_bin = $self->param_required("trimmomatic_bin");
  my $adapters = $self->param_required("trim_adapters_se");
  my $nthreads = $self->param("threads");

  my $cmd = "java -jar $trim_bin SE";
  $cmd .= " -threads $nthreads";
  $cmd .= " $seq";
  $cmd .= " $seq_trim";
  $cmd .= " ILLUMINACLIP:$adapters:2:30:10";
  $cmd .= " LEADING:3";
  $cmd .= " TRAILING:3";
  $cmd .= " SLIDINGWINDOW:4:15";
  $cmd .= " MINLEN:20";

  my ($stdout, $stderr, $exit) = capture {
    system($cmd);
  };

  if ($exit != 0) {
    die("Trimming failed: $stderr");
  }

  return $cmd;
}

sub trim_reads_paired {
  my ($self, $seq1, $seq2, $seq_trim1, $seq_trim2) = @_;
  return if -s $seq_trim1 and -s $seq_trim2;

  my $trim_bin = $self->param_required("trimmomatic_bin");
  my $adapters = $self->param_required("trim_adapters_pe");
  my $nthreads = $self->param("threads");

  my $unpaired_1 = $seq1 . "_unpaired";
  my $unpaired_2 = $seq2 . "_unpaired";

  my $cmd = "java -jar $trim_bin PE";
  $cmd .= " -threads $nthreads";
  $cmd .= " $seq1";
  $cmd .= " $seq2";
  $cmd .= " $seq_trim1";
  $cmd .= " $unpaired_1";
  $cmd .= " $seq_trim2";
  $cmd .= " $unpaired_2";
  $cmd .= " ILLUMINACLIP:$adapters:2:30:10";
  $cmd .= " LEADING:3";
  $cmd .= " TRAILING:3";
  $cmd .= " SLIDINGWINDOW:4:15";
  $cmd .= " MINLEN:20";

  my ($stdout, $stderr, $exit) = capture {
    system($cmd);
  };

  unlink $unpaired_1;
  unlink $unpaired_2;

  if ($exit != 0) {
    die("Trimming failed: $stderr");
  }

  return $cmd;
}

1;
