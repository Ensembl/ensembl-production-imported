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

package Bio::EnsEMBL::EGPipeline::Xref::DeleteUnattachedXref;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    'db_type' => 'core',
  };
}

sub run {
  my ($self) = @_;
  my $db_type = $self->param_required('db_type');

  my $dba = $self->get_DBAdaptor($db_type);
  my $dbh = $dba->dbc->db_handle();

  $self->drop_tmp_table($dbh);
  $self->create_tmp_table($dbh);
  $self->delete_xrefs($dbh);
  $self->drop_tmp_table($dbh);

  $dba->dbc and $dba->dbc->disconnect_if_idle();
}

sub drop_tmp_table {
  my ($self, $dbh) = @_;

  my $sql = "DROP TEMPORARY TABLE IF EXISTS xref_used";
  my $sth = $dbh->prepare($sql);
  $sth->execute();
}

sub create_tmp_table {
  my ($self, $dbh) = @_;

  my $sql = "
    CREATE TEMPORARY TABLE xref_used (INDEX (xref_id)) AS
    SELECT 'associated_xref' as table_name, xref_id
    FROM associated_xref
    INNER JOIN xref using (xref_id) GROUP BY xref_id
    UNION
    SELECT 'associated_xref' as table_name, source_xref_id
    FROM associated_xref
    INNER JOIN xref on source_xref_id = xref.xref_id GROUP BY source_xref_id
    UNION
    SELECT 'dependent_xref' as table_name, master_xref_id
    FROM dependent_xref
    INNER JOIN xref on master_xref_id = xref.xref_id GROUP BY master_xref_id
    UNION
    SELECT 'dependent_xref' as table_name, dependent_xref_id
    FROM dependent_xref
    INNER JOIN xref on dependent_xref_id = xref.xref_id GROUP BY dependent_xref_id
    UNION
    SELECT 'external_synonym' as table_name, xref_id
    FROM external_synonym
    INNER JOIN xref using (xref_id) GROUP BY xref_id
    UNION
    SELECT 'object_xref' as table_name, xref_id
    FROM object_xref
    INNER JOIN xref using (xref_id) GROUP BY xref_id
    UNION
    SELECT 'ontology_xref' as table_name, source_xref_id
    FROM ontology_xref
    INNER JOIN xref on source_xref_id = xref.xref_id GROUP BY source_xref_id
    UNION
    SELECT 'gene' as table_name, display_xref_id
    FROM gene
    INNER JOIN xref on display_xref_id = xref.xref_id GROUP BY display_xref_id
    UNION
    SELECT 'transcript' as table_name, display_xref_id
    FROM transcript
    INNER JOIN xref on display_xref_id = xref.xref_id GROUP BY display_xref_id
    UNION
    SELECT 'interpro' as table_name, xref_id
    FROM external_db
    INNER JOIN xref using (external_db_id) WHERE db_name = 'InterPro'
    GROUP BY xref_id
  ";

  my $sth = $dbh->prepare($sql);
  $sth->execute();
}

sub delete_xrefs {
  my ($self, $dbh) = @_;

  my $sql = "
    DELETE xref.* FROM xref LEFT OUTER JOIN xref_used USING (xref_id)
    WHERE xref_used.xref_id IS NULL
  ";
  my $sth = $dbh->prepare($sql);
  $sth->execute();
}

1;
