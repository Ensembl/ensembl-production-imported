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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::HtseqFactory;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Path::Tiny qw(path);

sub write_output {
  my ($self) = @_;

  my $bam = $self->param_required('bam_file');
  my $is_stranded = $self->param_required('is_stranded');

  my @strands = qw(firststrand secondstrand);
  my @numbers = qw(unique total);
  my @cases;
  # We need htseq-count for:
  if ($is_stranded) {
    for my $number (@numbers) {
      for my $strand (@strands) {
        push @cases, { bam_file => $bam, strand => $strand, number => $number };
      }
    }
  } else {
    for my $number (@numbers) {
      push @cases, { bam_file => $bam, strand => 'unstranded', number => $number };
    }
  }

  for my $case (@cases) {
    $self->dataflow_output_id($case, 2);
  }
}

1;
