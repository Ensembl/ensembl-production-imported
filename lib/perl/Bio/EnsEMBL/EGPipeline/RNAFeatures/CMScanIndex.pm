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

package Bio::EnsEMBL::EGPipeline::RNAFeatures::CMScanIndex;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable::CMScan;

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub run {
  my ($self) = @_;
  my $species         = $self->param_required('species');
  my $rfam_cm_file    = $self->param_required('rfam_cm_file');
  my $rfam_logic_name = $self->param_required('rfam_logic_name');
  my $cm_file         = $self->param_required('cmscan_cm_file');
  my $logic_name      = $self->param_required('cmscan_logic_name');
  my $parameters_hash = $self->param_required('parameters_hash');
  
  my $aa = $self->core_dba->get_adaptor('Analysis');
  
  if (! exists $$cm_file{$species} && ! exists $$cm_file{'all'}) {
    $$parameters_hash{'-cm_file'}  = $rfam_cm_file;
    $$parameters_hash{'-analysis'} = $aa->fetch_by_logic_name($rfam_logic_name);
  } else {
    $$parameters_hash{'-cm_file'} = $$cm_file{$species} || $$cm_file{'all'};
    $$parameters_hash{'-analysis'} = $aa->fetch_by_logic_name($logic_name);
  }
  
  my %parameters = %{$self->param('parameters_hash')};
  
  my $cmscan = Bio::EnsEMBL::Analysis::Runnable::CMScan->new(%parameters);
  
  $cmscan->prepare_index();

  my $db_connection = $aa->dbc();
  $db_connection && $db_connection->disconnect_if_idle();
}

1;
