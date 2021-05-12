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

Bio::EnsEMBL::EGPipeline::Xref::LoadUniProt

=head1 DESCRIPTION

Add UniProt xrefs to a core database, based on UniParc matches.

=head1 Author

Dan Staines and James Allen

=cut

package Bio::EnsEMBL::EGPipeline::Xref::LoadUniProt;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::Xref::LoadXref');

use Bio::EnsEMBL::DBSQL::DBAdaptor;

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    'logic_name'             => 'xrefuniparc',
    'external_dbs'           => {
                                  reviewed   => 'Uniprot/SWISSPROT',
                                  unreviewed => 'Uniprot/SPTREMBL',
                                  splicevar  => 'Uniprot/Varsplic',
                                },
    'replace_all'            => 0,
    'gene_name_source'       => [],
    'overwrite_gene_name'    => 0,
    'description_source'     => [],
    'overwrite_description'  => 0,
    'description_blacklist'  => [],
    'uniparc_external_db'    => 'UniParc',
    'uniprot_gn_external_db' => 'Uniprot_gn',
  };
}

sub run {
  my ($self) = @_;
  my $db_type      = $self->param_required('db_type');
  my $logic_name   = $self->param_required('logic_name');
  my $external_dbs = $self->param_required('external_dbs');
  my $replace_all  = $self->param_required('replace_all');

  my @external_dbs = values %$external_dbs;

  my $dba = $self->get_DBAdaptor($db_type);
  my $aa  = $dba->get_adaptor('Analysis');

  my $analysis = $aa->fetch_by_logic_name($logic_name);

  if ($replace_all) {
    $self->remove_xrefs($dba, \@external_dbs);
  }

  foreach my $external_db (@external_dbs) {
    $self->external_db_reset($dba, $external_db);
  }

  my $all_uniprots = $self->add_xrefs($dba, $analysis, $external_dbs);

  foreach my $external_db (@external_dbs) {
    $self->external_db_update($dba, $external_db);
  }

  $self->set_descriptions($dba, $analysis, $external_dbs);

  $self->set_gene_names($dba, $analysis, $all_uniprots);

  $dba->dbc and $dba->dbc->disconnect_if_idle();
}

sub add_xrefs {
  my ($self, $dba, $analysis, $external_dbs) = @_;
  my $uniparc_db          = $self->param_required('uniparc_db');
  my $uniprot_db          = $self->param_required('uniprot_db');
  my $uniparc_external_db = $self->param_required('uniparc_external_db');

  my $uniparc_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%$uniparc_db);
  my $uniprot_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%$uniprot_db);

  my @all_uniprots;

  my $ta   = $dba->get_adaptor('Translation');
  my $dbea = $dba->get_adaptor('DBEntry');

  my $tax_id = $dba->get_MetaContainer->get_taxonomy_id();

  my $translations = $ta->fetch_all();
  foreach my $translation (@$translations) {
    my $upi_xrefs = $dbea->fetch_all_by_Translation($translation, $uniparc_external_db);

    foreach my $upi_xref (@$upi_xrefs) {
      my $uniprots = $self->get_uniprot_for_upi($uniparc_dba, $uniprot_dba, $tax_id, $upi_xref->primary_id);

      foreach my $uniprot (@$uniprots) {
        if (defined $$uniprot{description} && $$uniprot{description} eq $translation->stable_id) {
          delete $$uniprot{description};
        }

        my $xref = $self->add_xref($uniprot, $analysis, $external_dbs);
        if (defined $xref) {
              $dbea->store($xref, $translation->dbID(), 'Translation', undef, $upi_xref);
        }
      }

      push @all_uniprots, @$uniprots;
    }
  }

  return \@all_uniprots;
}

sub add_xref {
  my ($self, $uniprot, $analysis, $external_dbs) = @_;

  my $xref;

  if (exists $$external_dbs{$$uniprot{type}}) {
    my $external_db = $$external_dbs{$$uniprot{type}};

    $xref = Bio::EnsEMBL::DBEntry->new(
      -PRIMARY_ID  => $$uniprot{ac},
      -DISPLAY_ID  => $$uniprot{ac},
      -DESCRIPTION => $$uniprot{description},
      -VERSION     => $$uniprot{version},
      -DBNAME      => $external_db,
      -INFO_TYPE   => 'DEPENDENT',
    );
       $xref->analysis($analysis);
  }

  return $xref;
}

sub get_uniprot_for_upi {
  my ($self, $uniparc_dba, $uniprot_dba, $tax_id, $upi) = @_;
  
  my @blacklist = @{$self->param_required('description_blacklist')};
  my $blacklist = join('|', @blacklist);

  my %uniprots;

  my $uniparc_sql = q/
    SELECT ac FROM
      uniparc.xref x INNER JOIN
      uniparc.cv_database d ON (x.dbid = d.id)
    WHERE
      taxid = ? AND
      upi = ? AND
      descr IN ('TREMBL','SWISSPROT', 'SWISSPROT_VARSPLIC') AND
      deleted='N'
  /;
  my $uniparc_sth = $uniparc_dba->dbc->db_handle->prepare($uniparc_sql);

  my $uniprot_sql = q/
    SELECT
      d.name,
      REPLACE(dd1.descr, '^') AS description,
      d.entry_type,
      s.version,
      gn.name AS gene_name,
      gnt.type AS gene_name_type
    FROM
      sptr.dbentry d JOIN
      sequence s ON (s.dbentry_id = d.dbentry_id) LEFT OUTER JOIN
      sptr.dbentry_2_desc dd1 ON (dd1.dbentry_id = d.dbentry_id) LEFT OUTER JOIN
      cv_desc cd1 ON (dd1.desc_id = cd1.desc_id) LEFT OUTER JOIN
      gene g ON (d.dbentry_id = g.dbentry_id) LEFT OUTER JOIN
      gene_name gn ON (gn.GENE_ID = g.GENE_ID) LEFT OUTER JOIN
      cv_gene_name_type gnt ON (gnt.GENE_NAME_TYPE_ID = gn.GENE_NAME_TYPE_ID)
    WHERE
      d.accession = ? AND
      (
        cd1.desc_id IS NULL OR
        (
          cd1.section_type = 'Main' AND
          cd1.catg_type IN ('RecName', 'SubName') AND
          cd1.subcatg_type = 'Full'
        )
      )
  /;
  my $uniprot_sth = $uniprot_dba->dbc->db_handle->prepare($uniprot_sql);

  $uniparc_sth->execute($tax_id, $upi);

  my $uniparc_results = $uniparc_sth->fetchall_arrayref();
  foreach my $uniparc_result (@$uniparc_results) {
    my $ac = $$uniparc_result[0];
    my %synonyms;

    $uniprot_sth->execute($ac);
    my $uniprot_results = $uniprot_sth->fetchall_arrayref();
    foreach my $uniprot_result (@$uniprot_results) {
      my ($name, $desc, $type_code, $version, $gene_name, $gene_name_type) = @$uniprot_result;
      
      if (!exists $uniprots{$ac}) {
        my $type;
        if ($type_code == 0) {
          $type = 'reviewed';
        } elsif ($type_code == 1) {
          $type = 'unreviewed';
        } elsif ($type_code == 13) {
          $type = 'splicevar';
        } else {
          $self->throw("Unrecognised UniProt type code '$type_code' for $ac");
        } 
        
        $uniprots{$ac}{ac}   = $ac;
        $uniprots{$ac}{name} = $ac;
        $uniprots{$ac}{type} = $type;
        if (defined $desc && $desc ne '') {
          if ($desc !~ /^($blacklist)$/) {
            $uniprots{$ac}{description} = $desc;
          }
        }
        if (defined $version && $version ne '') {
          $uniprots{$ac}{version} = $version;
        }
      }


      if (defined $gene_name) {
        if (!exists $synonyms{$gene_name}) {
          if ($gene_name_type eq 'Name' && $gene_name !~ /^\d+$/) {
            $uniprots{$ac}{gene_name} = $gene_name;
          } else {
            push @{$uniprots{$ac}{synonyms}}, $gene_name;
          }
          $synonyms{$gene_name}++;
        }
      }
    }
  }

  my @uniprots = values %uniprots;

  return \@uniprots;
}

sub set_descriptions {
  my ($self, $dba, $analysis, $external_dbs) = @_;
  my $sources   = $self->param_required('description_source');
  my $overwrite = $self->param_required('overwrite_description');
  
  if (@$sources) {
    my @db_names = map { $$external_dbs{$_} } @$sources;
    my $db_names = "'" . join("','", @db_names) . "'";

    foreach my $external_db (@db_names) {
      $self->remove_descriptions($dba, $external_db);
    }

    my $sql = q/
      UPDATE
        gene g INNER JOIN
        translation t ON g.canonical_transcript_id = t.transcript_id INNER JOIN
        object_xref ox ON (ox.ensembl_id = t.translation_id) INNER JOIN
        xref x USING (xref_id) INNER JOIN
        external_db edb USING (external_db_id)
      SET
        g.description = CONCAT(x.description, " [Source:", edb.db_display_name, ";Acc:", x.dbprimary_acc, "]")
      WHERE
        x.description IS NOT NULL AND
        ox.ensembl_object_type = 'Translation' AND
        ox.analysis_id = ? AND
    /;
    $sql .= " edb.db_name IN ($db_names) ";
    $sql .= " AND g.description IS NULL " unless $overwrite;

    my $sth = $dba->dbc->db_handle->prepare($sql);
    $sth->execute($analysis->dbID);
  }
}

sub remove_descriptions {
  my ($self, $dba, $external_db) = @_;

  my $sql = q/
    UPDATE
      gene g,
      external_db edb
    SET
      g.description = NULL
    WHERE
      edb.db_name = ? AND
      g.description LIKE CONCAT('%', ' [Source:', edb.db_display_name, '%')
  /;
  my $sth = $dba->dbc->db_handle->prepare($sql);
  $sth->execute($external_db);
}

sub set_gene_names {
  my ($self, $dba, $analysis, $all_uniprots) = @_;
  my $sources     = $self->param_required('gene_name_source');
  my $overwrite   = $self->param_required('overwrite_gene_name');
  my $external_db = $self->param_required('uniprot_gn_external_db');

  if (@$sources) {
    my $synonym_sql = 'INSERT IGNORE INTO external_synonym VALUES (?, ?)';
    my $synonym_sth = $dba->dbc->db_handle->prepare($synonym_sql);
  
    $self->remove_gene_names($dba, $external_db);

    my $dbea = $dba->get_adaptor('DBEntry');

    my %sources = map { $_ => 1 } @$sources;

    foreach my $uniprot (@$all_uniprots) {
      if (defined $$uniprot{gene_name} && exists $sources{$$uniprot{type}}) {
        my $xref = Bio::EnsEMBL::DBEntry->new(
          -PRIMARY_ID  => $$uniprot{gene_name},
          -DISPLAY_ID  => $$uniprot{gene_name},
          -DBNAME      => $external_db,
          -INFO_TYPE   => 'DEPENDENT',
        );

        $dbea->store($xref);
        
        # Synonyms are only added by the '$dbea->store' method if
        # the xref doesn't already exist...
        if (defined $$uniprot{synonyms}) {
          for my $synonym (@{$$uniprot{synonyms}}) {
            $synonym_sth->execute($xref->dbID, $synonym)
              or $self->throw("Failed to add synonym '$synonym'");
          }
        }
        
        my $sql = q/
          UPDATE
            gene g INNER JOIN
            translation t ON g.canonical_transcript_id = t.transcript_id INNER JOIN
            object_xref ox ON (ox.ensembl_id = t.translation_id) INNER JOIN
            xref x USING (xref_id) INNER JOIN
            external_db edb USING (external_db_id)
          SET
            g.display_xref_id = ?
          WHERE
            ox.ensembl_object_type = 'Translation' AND
            ox.analysis_id = ? AND
            x.dbprimary_acc = ?
        /;
        $sql .= " AND g.display_xref_id IS NULL " unless $overwrite;

        my $sth = $dba->dbc->db_handle->prepare($sql);
        $sth->execute($xref->dbID, $analysis->dbID, $uniprot->{ac});
      }
    }
  }
}

sub remove_gene_names {
  my ($self, $dba, $external_db) = @_;

  my $sql = q/
    UPDATE
      gene g INNER JOIN
      xref x ON g.display_xref_id = x.xref_id INNER JOIN
      external_db edb USING (external_db_id)
    SET
      g.display_xref_id = NULL
    WHERE
      edb.db_name = ?
  /;
  my $sth = $dba->dbc->db_handle->prepare($sql);
  $sth->execute($external_db);
}

1;
