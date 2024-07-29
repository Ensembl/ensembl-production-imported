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

sub param_defaults {
  my ($self) = @_;
  
  return {
    # Features on which to count the reads
    features => ['exon'],
  };
}

sub write_output {
  my ($self) = @_;

  my $bam = $self->param_required('bam_file');
  my $aligner_metadata = $self->param_required('aligner_metadata');
  my $is_stranded = $aligner_metadata->{'is_stranded'};
  my $features = $self->param_required('features');

  my @strands = qw(firststrand secondstrand);
  my @numbers = qw(unique total);
  my @cases;
  
  for my $feature (@$features) {
    # We need htseq-count for:
    if ($is_stranded) {
      for my $number (@numbers) {
        for my $strand (@strands) {
          push @cases, {
            bam_file => $bam,
            strand => $strand,
            number => $number,
            feature => $feature,
          };
        }
      }
    } else {
      for my $number (@numbers) {
        push @cases, {
          bam_file => $bam,
          strand => 'unstranded',
          number => $number,
          feature => $feature,
        };
      }
    }
  }

  for my $case (@cases) {
    $self->dataflow_output_id($case, 2);
  }
}

1;
