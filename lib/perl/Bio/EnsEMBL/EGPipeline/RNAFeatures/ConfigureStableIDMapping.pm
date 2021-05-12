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

Bio::EnsEMBL::EGPipeline::RNAFeatures::ConfigureStableIDMapping

=cut

package Bio::EnsEMBL::EGPipeline::RNAFeatures::ConfigureStableIDMapping;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub write_output {
  my ($self) = @_;
  my $species         = $self->param_required('species');
  my $all_new_species = $self->param_required('all_new_species');
  
  if (!$all_new_species) {
    my $old_reg_conf = $self->param_required('old_reg_conf');
    
    my $old_reg = 'Bio::EnsEMBL::Registry';
    $old_reg->load_all($old_reg_conf);
    my $old_dba = $old_reg->get_DBAdaptor($species, 'core');
    
    if (defined $old_dba) {
      my $old_db_params = {
        old_host   => $old_dba->dbc->host,
        old_port   => $old_dba->dbc->port,
        old_user   => $old_dba->dbc->user,
        old_pass   => $old_dba->dbc->pass,
        old_dbname => $old_dba->dbc->dbname,
      };
      
      $self->dataflow_output_id($old_db_params, 2);
    } else {
      $self->warning("No old database found for $species, therefore no stable ID mapping was done.");
    }
  }
}

1;
