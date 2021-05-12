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

package Bio::EnsEMBL::EGPipeline::Xref::ConfigureSource;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub write_output {
  my ($self) = @_;
  my $logic_name  = $self->param_required('logic_name');
  my $external_db = $self->param_required('external_db');

  my ($database_type, $query_type, $query_subtype, $data_source);
  
  if ($external_db =~ /^Uniprot/) {
    ($database_type, $query_type, $query_subtype, $data_source) = ('pep', 'pep', '', 'uniprot');
  } elsif ($external_db eq 'RefSeq_peptide') {
    ($database_type, $query_type, $query_subtype, $data_source) = ('pep', 'pep', '', 'refseq');
  } elsif ($external_db eq 'RefSeq_dna') {
    ($database_type, $query_type, $query_subtype, $data_source) = ('dna', 'dna', 'transcript', 'refseq');
  }
  
  $self->dataflow_output_id(
    {
      'logic_name'    => $logic_name,
      'external_db'   => $external_db,
      'database_type' => $database_type,
      'query_type'    => $query_type,
      'query_subtype' => $query_subtype,
      'data_source'   => $data_source,
    }, 1);
}

1;
