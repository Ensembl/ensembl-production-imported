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

Bio::EnsEMBL::EGPipeline::Xref::LoadUniProtGO

=head1 DESCRIPTION

Add UniProt transitive GO xrefs to a core database.

=head1 Author

Dan Staines and James Allen

=cut

package Bio::EnsEMBL::EGPipeline::Xref::LoadUniProtGO;

use strict;
use warnings;
use feature 'say';

use base ('Bio::EnsEMBL::EGPipeline::Xref::LoadXref');

use Bio::EnsEMBL::DBSQL::DBAdaptor;

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    'logic_name'           => 'gouniprot',
    'external_db'          => 'GO',
    'uniprot_external_dbs' => {reviewed => 'Uniprot/SWISSPROT', unreviewed => 'Uniprot/SPTREMBL'},
    'replace_all'          => 0,
  };
}

sub run {
  my ($self) = @_;
  my $db_type     = $self->param_required('db_type');
  my $logic_name  = $self->param_required('logic_name');
  my $external_db = $self->param_required('external_db');
  my $replace_all = $self->param_required('replace_all');

  my $dba = $self->get_DBAdaptor($db_type);
  my $aa  = $dba->get_adaptor('Analysis');

  my $analysis = $aa->fetch_by_logic_name($logic_name);

  if ($replace_all) {
    $self->remove_xrefs($dba, [$external_db]);
  }

  $self->external_db_reset($dba, $external_db);

  $self->add_xrefs($dba, $analysis, $external_db);

  $self->external_db_update($dba, $external_db);

  $dba->dbc and $dba->dbc->disconnect_if_idle();
}

sub add_xrefs {
  my ($self, $dba, $analysis, $external_db) = @_;
  my $uniprot_db           = $self->param_required('uniprot_db');
  my $uniprot_external_dbs = $self->param_required('uniprot_external_dbs');

  my $uniprot_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%$uniprot_db);

  my $ta   = $dba->get_adaptor('Translation');
  my $dbea = $dba->get_adaptor('DBEntry');

  my $translations = $ta->fetch_all();
  foreach my $translation (@$translations) {
    foreach my $status (keys %$uniprot_external_dbs) {
      my $uniprot_xrefs = $dbea->fetch_all_by_Translation($translation, $$uniprot_external_dbs{$status});

      foreach my $uniprot_xref (@$uniprot_xrefs) {
        if ($uniprot_xref->info_type eq 'DEPENDENT') {
          my $gos = $self->get_go_for_uniprot($uniprot_dba, $uniprot_xref->primary_id);

          foreach my $go (@$gos) {
            my $xref = $self->add_xref($go, $uniprot_xref, $analysis, $external_db);
            $dbea->store($xref, $translation->transcript->dbID(), 'Transcript');
          }
        }
      }
    }
  }
}

sub add_xref {
  my ($self, $go, $uniprot_xref, $analysis, $external_db) = @_;

  my ($term, $description, $evidence) = @$go;

  my $xref = Bio::EnsEMBL::OntologyXref->new(
    -PRIMARY_ID  => $term,
    -DISPLAY_ID  => $term,
    -DESCRIPTION => $description,
    -DBNAME      => $external_db,
    -INFO_TYPE   => 'DEPENDENT',
  );
  $xref->analysis($analysis);
  $xref->{version} = undef;

  if ($evidence) {
    $xref->add_linkage_type($evidence, $uniprot_xref);
  }

  return $xref;
}

sub get_go_for_uniprot {
  my ($self, $uniprot_dba, $accession) = @_;

  my $sql = q/
    SELECT
      primary_id AS term,
      regexp_replace(secondary_id,'.*:','') AS description,
      regexp_replace(note,':.*','') AS evidence
    FROM
      dbentry d INNER JOIN
      dbentry_2_database dd USING (dbentry_id)
    WHERE
      d.accession = ? AND
      dd.database_id = 'GO'
  /;
  my $sth = $uniprot_dba->dbc->db_handle->prepare($sql);

  $sth->execute($accession);

  my $gos = $sth->fetchall_arrayref();

  return $gos;
}

1;
