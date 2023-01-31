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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::MetadataCheck;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub run {
  my ($self) = @_;
  my $seq1 = $self->param_required('seq_file_1');
  my $seq2 = $self->param('seq_file_2');
  my $input_metadata = $self->param('input_metadata');

  if (defined $seq1 = $self->param_required('seq_file_1') ){
    print $input_metadata;
  }
    #$output->{subset_seq_file_2} = $self->param('subset_seq_file_2');
 #} else {
  #  $output->{subset_seq_file_2} = undef;

  # Check that there are genes!
  #my $merged_data =
 # if count ($merged_data) > 1{
   # paired
  }
1;