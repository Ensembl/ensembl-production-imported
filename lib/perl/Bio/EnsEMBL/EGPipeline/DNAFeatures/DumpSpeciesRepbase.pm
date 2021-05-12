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

package Bio::EnsEMBL::EGPipeline::DNAFeatures::DumpSpeciesRepbase;

use strict;
use warnings;

use Bio::EnsEMBL::EGPipeline::DNAFeatures::Utils qw(get_repbase_species_lib);

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub run {
  my ($self) = @_;
  
  my $species_repeat_lib = $self->param_required('species_repeat_lib');
  
  my $dba = $self->get_DBAdaptor('core');
  my $species = get_repbase_species_lib($dba, $species_repeat_lib);
}

1;
