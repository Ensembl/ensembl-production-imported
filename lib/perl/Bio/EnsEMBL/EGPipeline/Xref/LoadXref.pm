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

Bio::EnsEMBL::EGPipeline::Xref::LoadXref

=head1 DESCRIPTION

Generic class for adding xrefs to a core database.

=head1 Author

James Allen

=cut

use strict;
use warnings;

package Bio::EnsEMBL::EGPipeline::Xref::LoadXref;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');
use Time::Piece;

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    'db_type' => 'core',
  };
}

sub external_db_reset {
  my ($self, $dba, $db_name) = @_;

  my $dbh = $dba->dbc->db_handle();
  my $sql = "UPDATE external_db SET db_release = NULL WHERE db_name = ?;";
  my $sth = $dbh->prepare($sql);
  $sth->execute($db_name) or $self->throw("Failed to execute ($db_name): $sql");
}

sub external_db_update {
  my ($self, $dba, $db_name) = @_;

  my $t = localtime;
  my $db_release = "EG Xref pipeline; ".$t->datetime;

  my $dbh = $dba->dbc->db_handle();
  my $sql = "UPDATE external_db SET db_release = ? WHERE db_name = ?;";
  my $sth = $dbh->prepare($sql);
  $sth->execute($db_release, $db_name) or $self->throw("Failed to execute ($db_release, $db_name): $sql");
}

sub remove_xrefs {
  my ($self, $dba, $external_dbs) = @_;

  my $db_names = "'" . join("','", @$external_dbs) . "'";

  my $sql = q/
    DELETE ox.*, ix.*, ontix.* FROM
      translation tn INNER JOIN
      transcript tt USING (transcript_id) INNER JOIN
      seq_region sr USING (seq_region_id) INNER JOIN
      coord_system cs USING (coord_system_id) INNER JOIN
      object_xref ox ON (ox.ensembl_id = tn.translation_id) INNER JOIN
      xref x USING (xref_id) INNER JOIN
      external_db edb USING (external_db_id) LEFT OUTER JOIN
      identity_xref ix USING (object_xref_id) LEFT OUTER JOIN
      ontology_xref ontix USING (object_xref_id)
    WHERE
      ox.ensembl_object_type = 'Translation' AND
      cs.species_id = ? AND
  /;
  $sql .= "edb.db_name IN ($db_names)";
  my $sth = $dba->dbc->db_handle->prepare($sql);
  $sth->execute($dba->species_id);

  $sql = q/
    DELETE ox.*, ix.*, ontix.* FROM
      transcript tt INNER JOIN
      seq_region sr USING (seq_region_id) INNER JOIN
      coord_system cs USING (coord_system_id) INNER JOIN
      object_xref ox ON (ox.ensembl_id = tt.transcript_id) INNER JOIN
      xref x USING (xref_id) INNER JOIN
      external_db edb USING (external_db_id) LEFT OUTER JOIN
      identity_xref ix USING (object_xref_id) LEFT OUTER JOIN
      ontology_xref ontix USING (object_xref_id)
    WHERE
      ox.ensembl_object_type = 'Transcript' AND
      cs.species_id = ? AND
  /;
  $sql .= "edb.db_name IN ($db_names)";
  $sth = $dba->dbc->db_handle->prepare($sql);
  $sth->execute($dba->species_id);

  $sql = q/
    DELETE ox.*, dx.*, ix.*, ontix.* FROM
      gene g INNER JOIN
      seq_region sr USING (seq_region_id) INNER JOIN
      coord_system cs USING (coord_system_id) INNER JOIN
      object_xref ox ON (ox.ensembl_id = g.gene_id) INNER JOIN
      xref x USING (xref_id) INNER JOIN
      external_db edb USING (external_db_id) LEFT OUTER JOIN
      dependent_xref dx USING (object_xref_id) LEFT OUTER JOIN
      identity_xref ix USING (object_xref_id) LEFT OUTER JOIN
      ontology_xref ontix USING (object_xref_id)
    WHERE
      ox.ensembl_object_type = 'Gene' AND
      cs.species_id = ? AND
  /;
  $sql .= "edb.db_name IN ($db_names)";
  $sth = $dba->dbc->db_handle->prepare($sql);
  $sth->execute($dba->species_id);
}

1;
