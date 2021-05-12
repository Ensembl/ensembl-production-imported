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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::DatasetSpeciesFactory;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use JSON;

sub run {
  my ($self) = @_;
  
  my $species = $self->param("species");
  
  my $dba = $self->core_dba();
  my $ma = $dba->get_adaptor('MetaContainer');
  my ($organism) = @{ $ma->list_value_by_key("BRC4.organism_abbrev") };
  my ($component) = @{ $ma->list_value_by_key("BRC4.component") };
  
  # Slurp json
  my $json;
  {
    local $/; #enable slurp
    open my $fh, "<", $self->param_required('datasets_file');
    $json = <$fh>;
  }
  my $data = decode_json($json);

  # Check there is at least one dataset for this species
  my ($dataset) = grep { $_->{species} eq $organism } @$data;
  if ($dataset) {
    my $data = {
      species => $species,
      component => $component,
      organism => $organism
    };
    $self->dataflow_output_id($data, 2);
  }
}

1;
