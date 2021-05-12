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

    Bio::EnsEMBL::EGPipeline::PHIBase::XrefsCleaner

=head1 SYNOPSIS

   
    
=head1 DESCRIPTION

    This is a generic RunnableDB module for cleaning all PHI-base associated xrefs from all the dbs in the registry for a particular division.
    It's a copy/paste implementation of ensembl-production's module Bio::EnsEMBL::Production::Pipeline::Common::DbCmd.pm which
    executes a set of SQL commands from a given file by running mysql on the command line.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <https://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <https://www.ensembl.org/Help/Contact>.

=cut


package Bio::EnsEMBL::EGPipeline::PHIBase::XrefsCleaner;

use strict;
use warnings;

use base (
  'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
  'Bio::EnsEMBL::Production::Pipeline::Common::Base',
  'Bio::EnsEMBL::Hive::Process'
);

sub param_defaults {
  my ($self) = @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'db_type' => 'core',
    'input_query'  => 'DELETE x.*,ox.*,ax.*,ag.*,oox.* 
                        from external_db e
                        join xref x using (external_db_id)
                        join object_xref ox using (xref_id)
                        join associated_xref ax using (object_xref_id)
                        join associated_group ag using (associated_group_id)
                        join ontology_xref oox using (object_xref_id)
                        WHERE e.db_name="PHI";

                        DELETE from gene_attrib WHERE attrib_type_id=358;

                        DELETE from gene_attrib WHERE attrib_type_id=317 and value="PHI";',
  };
}

sub fetch_input {
  my $self = shift @_;
  my $db_type = $self->param('db_type');
  
  if ($db_type eq 'hive') {
    $self->param('db_conn', $self->dbc);
  } else {
    $self->param('db_conn', $self->get_DBAdaptor($db_type)->dbc);
  }
  my $dbname = $self->param('dbname');
  my $species = $self->param('species');
  if ($dbname && $dbname ne '') {
    print "------------------------------\nDeleting existing phibase xrefs from $dbname \n";
  }

  $self->SUPER::fetch_input();
  
}

1;
