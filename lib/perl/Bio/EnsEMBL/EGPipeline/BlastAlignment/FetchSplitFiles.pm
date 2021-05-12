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

Bio::EnsEMBL::EGPipeline::BlastAlignment::FetchSplitFiles

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::BlastAlignment::FetchSplitFiles;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub write_output {
  my ($self) = @_;
  my $species  = $self->param_required('species');
  my $seq_type = $self->param_required('seq_type');
  
  my $sql = "SELECT split_file FROM split_$seq_type WHERE species = ?;";
  my $sth = $self->hive_dbh->prepare($sql);
  $sth->execute($species);
  
  my $results = $sth->fetchall_arrayref();
  foreach my $result (@$results) {
    $self->dataflow_output_id({'split_file' => $$result[0]}, 2);
  }
  
  $self->dataflow_output_id({'seq_type' => $seq_type}, 1);
}

1;
