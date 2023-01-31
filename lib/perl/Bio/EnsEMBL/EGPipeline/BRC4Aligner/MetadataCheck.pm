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
  my $input_is_paired = $input_metadata->{'is_paired'};
  my $input_is_stranded = $input_metadata->{'is_stranded'};

  
  if (defined $seq1 && $seq2){
    $input_is_paired = "1" ;
  }
  else{
   print "It is not paired";
  }

  my $alter_metadata = {
      is_stranded => $input_is_stranded,
      is_paired => $input_is_paired,
  };
  
  $self->dataflow_output_id({
      alter_metadata => $alter_metadata
    },  1);
  }
1;