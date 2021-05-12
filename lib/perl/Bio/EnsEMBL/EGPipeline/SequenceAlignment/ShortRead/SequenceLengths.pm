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

package Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::SequenceLengths;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::SeqIO;

sub run {
  my ($self) = @_;
  
  my $fasta_file  = $self->param_required('fasta_file');
  my $length_file = $self->param_required('length_file');
  
  open (my $fh, '>', $length_file) or die "Failed to open file '$length_file'";
  
  my $seqs = Bio::SeqIO->new(-format => 'Fasta', -file => $fasta_file);
  while (my $seq = $seqs->next_seq) {
    print $fh $seq->id."\t".$seq->length."\n";
  }
  
  close($fh);
}

1;
