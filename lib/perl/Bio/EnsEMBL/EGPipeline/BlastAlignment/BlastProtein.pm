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


=pod

=head1 NAME

Bio::EnsEMBL::EGPipeline::BlastAlignment::BlastProtein

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::BlastAlignment::BlastProtein;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub param_defaults {
  return {
    'proteome_source' => 'file',
  };
}

sub write_output {
  my ($self) = @_;
  my $proteome_source = $self->param_required('proteome_source');
  
  my ($output_ids, $flow);
  
  if ($proteome_source eq 'file') {
    $$output_ids{'db_fasta_file'} = $self->param_required('db_fasta_file');
    $flow = 2;
    
  } elsif ($proteome_source eq 'database') {
    $flow = 3;
    
  } elsif ($proteome_source eq 'uniprot') {
    $flow = 4;
    
  } else {
    $self->throw("Unrecognised proteome_source '$proteome_source'");
    
  }
  
  $self->dataflow_output_id($output_ids, $flow);
}

1;
