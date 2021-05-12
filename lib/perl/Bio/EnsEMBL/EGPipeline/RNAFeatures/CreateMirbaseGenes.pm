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

Bio::EnsEMBL::EGPipeline::RNAFeatures::CreateMirbaseGenes

=cut

package Bio::EnsEMBL::EGPipeline::RNAFeatures::CreateMirbaseGenes;

use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;

use base ('Bio::EnsEMBL::EGPipeline::RNAFeatures::CreateGenes');

sub run {
  my ($self) = @_;
  my $source_logic_name = $self->param_required('source_logic_name');
  my $id_db             = $self->param_required('id_db');
  
  my $dba  = $self->core_dba();
  my $dafa = $dba->get_adaptor('DnaAlignFeature');
  my $aa   = $dba->get_adaptor('Analysis');
  my $ga   = $dba->get_adaptor('Gene');
  my $dbea = $dba->get_adaptor('DBEntry');
  my $ida  = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%$id_db);
  
  my $gene_source = $self->generate_source($dba);
  
  my @features = @{ $dafa->fetch_all_by_logic_name($source_logic_name) };
  
  my $total_count = 0;
  
  FEATURE: foreach my $feature (@features) {
    if (! defined $feature->db_name) {
      my $dbid = $feature->dbID;
      $self->throw("Feature (ID=$dbid) lacks an external_db reference.");
    }
    
    $self->create_gene($dba, $aa, $ga, $dbea, $ida, $gene_source, $feature);
    $total_count++;
  }
  
  my $msg = "$total_count miRBase genes were created.";
  $self->warning($msg);
}

1;
