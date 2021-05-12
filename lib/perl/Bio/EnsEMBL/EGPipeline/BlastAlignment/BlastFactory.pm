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

package Bio::EnsEMBL::EGPipeline::BlastAlignment::BlastFactory;

use strict;
use warnings;

use Bio::SeqIO;

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub param_defaults {
  return {
    'max_seq_length' => 10000000,
  };
}

sub write_output {
  my ($self) = @_;
  
  my $species   = $self->param_required('species');
  my $queryfile = $self->param_required('queryfile');
  
  if (!-e $queryfile) {
    $self->throw("Query file '$queryfile' does not exist");
  } else {
    my $max_seq_length = $self->param('max_seq_length');
    if (defined $max_seq_length) {
      my $total_length = 0;
      my $fasta = Bio::SeqIO->new(-format => 'Fasta', -file => $queryfile);
      while (my $seq = $fasta->next_seq) {
        $total_length += $seq->length;
      }
      if ($total_length > $max_seq_length) {
        $fasta = Bio::SeqIO->new(-format => 'Fasta', -file => $queryfile);
        while (my $seq = $fasta->next_seq) {
          my $seq_length = $seq->length;
          my ($start, $end) = (1, $max_seq_length);
          $end = $seq_length if $end > $seq_length;
          
          while ($start <= $seq_length) {
            my $querylocation = $seq->id.":$start-$end";
            
            $self->dataflow_output_id(
              {
                'queryfile'     => undef,
                'querylocation' => $querylocation,
              }, 2);
          
            $start = $end + 1;
            $end = $start + $max_seq_length - 1;
            $end = $seq_length if $end > $seq_length;
          }
        }
      
      } else {
        $self->dataflow_output_id(
          {
            'queryfile' => $queryfile,
          }, 2);
      }
      
    } else {
      $self->dataflow_output_id(
        {
          'queryfile' => $queryfile,
        }, 2);
    }
  }
}

1;
