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
  
  # Get datasets data
  my $json;
  {
    local $/; #enable slurp
    open my $fh, "<", $self->param_required('datasets_file');
    $json = <$fh>;
  }
  my $data = decode_json($json);
  my %organisms = map { $_->{species} => 1 } @$data;
  my %prod_names = map { $_->{production_name} => 1 } @$data;
  
  # Get registry species
  my $reg = 'Bio::EnsEMBL::Registry';
  if ($self->param_is_defined('registry_file')) {
    $reg->load_all($self->param('registry_file'));
  }
  my $all_dbas = $reg->get_all_DBAdaptors(-GROUP => 'core');
  for my $dba (@$all_dbas) {
    my $ma = $dba->get_adaptor('MetaContainer');
    my ($organism) = @{ $ma->list_value_by_key("BRC4.organism_abbrev") };
    my ($component) = @{ $ma->list_value_by_key("BRC4.component") };
    my ($prod_name) = $dba->species;
    
    if (defined $organisms{$organism}) {
      my $data = {
        species => $prod_name,
        component => $component,
        organism => $organism
      };
      $self->dataflow_output_id($data, 2);
      delete $organisms{$organism};
    } elsif (defined $prod_names{$prod_name}) {
      my $data = {
        species => $prod_name,
        organism => $prod_name,
        component => 'component'
      };
      $self->dataflow_output_id($data, 2);
    }
    $dba->dbc->disconnect_if_idle();
  }
  
  # Flows the organisms not found in the registry to another place
  if (%organisms) {
    for my $organism (sort keys %organisms) {
      $self->dataflow_output_id({ organism => $organism }, 3);
    }
  }
}

1;
