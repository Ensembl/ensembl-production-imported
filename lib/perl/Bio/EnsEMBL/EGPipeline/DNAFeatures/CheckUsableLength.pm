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

package Bio::EnsEMBL::EGPipeline::DNAFeatures::CheckUsableLength;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub param_defaults {
  my ($self) = @_;
  
  return {
    min_slice => 5000,
    repeat_modeler_length => 40_000,
  };
}

sub run {
  my ($self) = @_;

  my $species = $self->param('species');
  
  # Flows to 2 if there is enough length to run the modeler
  if ($self->has_usable_length()) {
    $self->dataflow_output_id({}, 2);
  } else {
    $self->dataflow_output_id({ species => $species, reason => 'no_min_length' }, 3);
  }
}

# Repeat modeler requires a minimum total length to be run: check that we have enough in the database
sub has_usable_length {
  my ($self) = @_;

  my $min_slice = $self->param('min_slice_length');
  my $min_length = $self->param('repeat_modeler_length');
  
  my $dba = $self->get_DBAdaptor('core');
  my $sql = "SELECT sum(length) FROM seq_region WHERE length >= ?";
  
  my $dbh = $dba->dbc->db_handle;
  my $sth = $dbh->prepare($sql);
  $sth->execute($min_slice);
  
  my ($sum) = $sth->fetchrow_array;
  
  return $sum >= $min_length;
}

1;
