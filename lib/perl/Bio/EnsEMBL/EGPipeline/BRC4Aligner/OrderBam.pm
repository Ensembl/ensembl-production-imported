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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::OrderBam;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Capture::Tiny ':all';
use Path::Tiny qw(path);
use File::Basename qw(dirname);
use File::Spec::Functions qw(catdir);

sub write_output {
  my ($self) = @_;
  
  my $threads = $self->param("threads") // 1;
  my $memory = $self->param("memory") // 8000;

  my $input_bam = $self->param_required('input_bam');
  my $sorted_bam = $self->param_required('sorted_bam');
  my $aligner_metadata = $self->param_required('aligner_metadata');
  my $is_paired = $aligner_metadata->{'is_paired'};

  # Filter out unmapped reads or pairs
  my $filter = "";
  if (defined $is_paired) {
    if ($is_paired) {
      $filter = "-G 12";
    } else {
      $filter = "-F 4";
    }
  }
  my $convert_cmd = "samtools view -bS $input_bam $filter";

  # Second part: sorting the BAM
  # Calculate the correct memory per thread
  my $mem = $memory / $threads;
  
  # Samtools sort is too greedy: we give it less
  $mem *= 0.8;
  my $mem_limit = $mem . 'M';
  
  # Final sort command, BY NAME for HTSeq-count
  my $sort_cmd = "samtools sort -n -@ $threads -m $mem_limit -o $sorted_bam -O 'bam' -T $sorted_bam.sorting -";
  
  my $cmd = "$convert_cmd | $sort_cmd";
  
  # Run command
  my ($stdout, $stderr, $exit) = capture {
    system( $cmd );
  };
  die "Cannot execute $cmd: $stderr" if $exit != 0;
}

1;
