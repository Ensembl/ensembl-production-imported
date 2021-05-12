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

Bio::EnsEMBL::EGPipeline::Xref::LoadUniProtXrefs

=head1 DESCRIPTION

Add UniProt transitive xrefs to a core database.

=head1 Author

Dan Staines and James Allen

=cut

package Bio::EnsEMBL::EGPipeline::Xref::LoadUniProtXrefs;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::Xref::LoadXref');

use Bio::EnsEMBL::DBSQL::DBAdaptor;

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    'logic_name'             => 'xrefuniprot',
    'external_dbs'           => {},
    'secondary_external_dbs' => {},
    'uniprot_external_dbs'   => {reviewed => 'Uniprot/SWISSPROT', unreviewed => 'Uniprot/SPTREMBL'},
    'replace_all'            => 0,
  };
}

sub run {
  my ($self) = @_;
  my $db_type                = $self->param_required('db_type');
  my $logic_name             = $self->param_required('logic_name');
  my $external_dbs           = $self->param_required('external_dbs');
  my $secondary_external_dbs = $self->param_required('secondary_external_dbs');
  my $replace_all            = $self->param_required('replace_all');
  
  foreach my $secondary (values %$secondary_external_dbs) {
    $$external_dbs{$secondary} = $secondary;
  }
  my @external_dbs = values %$external_dbs;
  push @external_dbs, 'protein_id' if exists($$external_dbs{'EMBL'});
  
  my $dba = $self->get_DBAdaptor($db_type);
  my $aa  = $dba->get_adaptor('Analysis');

  my $analysis = $aa->fetch_by_logic_name($logic_name);

  if ($replace_all) {
    $self->remove_xrefs($dba, \@external_dbs);
  }

  foreach my $external_db (@external_dbs) {
    $self->external_db_reset($dba, $external_db);
  }

  $self->add_xrefs($dba, $analysis, $external_dbs);

  foreach my $external_db (@external_dbs) {
    $self->external_db_update($dba, $external_db);
  }

  $dba->dbc and $dba->dbc->disconnect_if_idle();
}

sub add_xrefs {
  my ($self, $dba, $analysis, $external_dbs) = @_;
  my $uniprot_db             = $self->param_required('uniprot_db');
  my $uniprot_external_dbs   = $self->param_required('uniprot_external_dbs');
  my $secondary_external_dbs = $self->param_required('secondary_external_dbs');

  my $uniprot_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%$uniprot_db);

  my $ta   = $dba->get_adaptor('Translation');
  my $dbea = $dba->get_adaptor('DBEntry');

  my $translations = $ta->fetch_all();
  foreach my $translation (@$translations) {
    foreach my $status (keys %$uniprot_external_dbs) {
      my $uniprot_xrefs = $dbea->fetch_all_by_Translation($translation, $$uniprot_external_dbs{$status});

      foreach my $uniprot_xref (@$uniprot_xrefs) {
        if ($uniprot_xref->info_type eq 'DEPENDENT') {
          my $transitive_xrefs = $self->get_xrefs_for_uniprot($uniprot_dba, $uniprot_xref->primary_id);

          foreach my $transitive_xref (@$transitive_xrefs) {
            my $xref = $self->add_xref($transitive_xref, $analysis, $external_dbs);
            if (defined $xref) {
              if ($xref->dbname eq 'EntrezGene' || $xref->dbname eq 'UniGene') {
                $dbea->store($xref, $translation->transcript->get_Gene->dbID(), 'Gene', 0, $uniprot_xref);
              } else {
                 $dbea->store($xref, $translation->dbID(), 'Translation', 0, $uniprot_xref);
                if (exists($$secondary_external_dbs{$xref->dbname})) {
                  if ($$transitive_xref[2] && $$transitive_xref[2] ne '-') {
                    $$transitive_xref[0] = $$secondary_external_dbs{$xref->dbname};
                    $$transitive_xref[1] = $$transitive_xref[2];
                    my $sec_xref = $self->add_xref($transitive_xref, $analysis, $external_dbs);
                    if ($sec_xref->dbname eq 'RefSeq_dna') {
                      $dbea->store($sec_xref, $translation->transcript->dbID(), 'Transcript', 0, $uniprot_xref);
                    } else {
                       $dbea->store($sec_xref, $translation->dbID(), 'Translation', 0, $uniprot_xref);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

sub add_xref {
  my ($self, $transitive_xref, $analysis, $external_dbs) = @_;

  my ($dbname, $primary_id, $secondary_id, $quaternary_id) = @$transitive_xref;
  my $xref;
  
  if (exists $$external_dbs{$dbname}) {
    my $external_db = $$external_dbs{$dbname};

    # For ENA we don't want genomic references where we have the CDS
    if ($dbname eq 'EMBL'    &&
      defined $secondary_id  && $secondary_id ne '-' &&
      defined $quaternary_id && $quaternary_id eq 'Genomic_DNA')
    {
      $primary_id  = $secondary_id;
      $external_db = 'protein_id';
    }

    $xref = Bio::EnsEMBL::DBEntry->new(
      -PRIMARY_ID  => $primary_id,
      -DISPLAY_ID  => $primary_id,
      -DBNAME      => $external_db,
      -INFO_TYPE   => 'DEPENDENT',
    );
    $xref->analysis($analysis);
  }

  return $xref;
}

sub get_xrefs_for_uniprot {
  my ($self, $uniprot_dba, $accession) = @_;

  my $sql = q/
    SELECT
      abbreviation AS dbname,
      primary_id,
      secondary_id,
      quaternary_id
    FROM
      dbentry d INNER JOIN
      dbentry_2_database dd USING (dbentry_id) INNER JOIN
      database_name db USING (database_id)
    WHERE
      d.accession = ?
  /;
  my $sth = $uniprot_dba->dbc->db_handle->prepare($sql);

  $sth->execute($accession);

  my $xrefs = $sth->fetchall_arrayref();

  return $xrefs;
}

1;
