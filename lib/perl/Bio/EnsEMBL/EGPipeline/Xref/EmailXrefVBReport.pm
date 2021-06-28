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

Bio::EnsEMBL::EGPipeline::Xref::EmailXrefVBReport

=head1 DESCRIPTION

Run a few useful queries on the xrefs, for a given species.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::Xref::EmailXrefVBReport;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EmailReport');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    'db_type' => 'core',
  };
}

sub fetch_input {
  my ($self) = @_;
  my $species              = $self->param_required('species');
  my $logic_name           = $self->param_required('logic_name');
  my $vb_external_db       = $self->param_required('vb_external_db');
  my $citation_external_db = $self->param_required('citation_external_db');
  
  my $dba = $self->get_DBAdaptor($self->param('db_type'));
  my $dbh = $dba->dbc->db_handle;

  my $reports = "The xref pipeline for $species has completed.\n";

  $reports .= $self->config_summary([$vb_external_db, $citation_external_db]);

  $reports .= "Summaries are below; note that the last two include pre-existing data.\n";

  $reports .= $self->xref_summary($dbh, $logic_name, 'Direct xrefs:');

  $reports .= $self->gene_data_summary($self->hive_dbh, $species, 'gene_descriptions', 'Gene description sources');
  $reports .= $self->gene_data_summary($self->hive_dbh, $species, 'gene_names', 'Gene name sources');

  $reports .= $self->xref_total_summary($dbh, 'All xrefs, pre-existing and newly-added:');
  $reports .= $self->xref_ontology_summary($dbh, 'Ontology xrefs, pre-existing and newly-added:');

  $self->param('text', $reports);
}

sub config_summary {
  my ($self, $xref_sources) = @_;

  my $summary = "Cross references have been loaded for the following databases:\n";
  $summary   .= join(', ', @$xref_sources)."\n";

  return $summary;
}

sub xref_summary {
  my ($self, $dbh, $logic_name, $title) = @_;

  my $sql = "
    SELECT
      db_name,
      COUNT(DISTINCT xref_id) AS xref_count,
      ensembl_object_type,
      COUNT(DISTINCT ensembl_id) AS ensembl_object_count
    FROM
      analysis INNER JOIN
      object_xref USING (analysis_id) INNER JOIN
      xref USING (xref_id) INNER JOIN
      external_db USING (external_db_id)
    WHERE
      logic_name = ?
    GROUP BY
      db_name,
      ensembl_object_type
    ;
  ";

  my $sth = $dbh->prepare($sql);
  $sth->execute($logic_name);

  my $columns = $sth->{NAME};
  my $results = $sth->fetchall_arrayref();

  return $self->format_table($title, $columns, $results);
}

sub gene_data_summary {
  my ($self, $dbh, $species, $table_name, $title) = @_;

  my $sql = "
    SELECT db_name, total FROM $table_name
    WHERE species = ? AND timing = ?
    ORDER BY db_name
    ;
  ";
  my $sth = $dbh->prepare($sql);

  $sth->execute($species, 'before');
  my $columns = $sth->{NAME};
  my $before_results = $sth->fetchall_arrayref();
  my $before_table = $self->format_table("Before pipeline:", $columns, $before_results);

  $sth->execute($species, 'after');
  my $after_results = $sth->fetchall_arrayref();
  my $after_table = $self->format_table("After pipeline:", $columns, $after_results);

  return "\n$title$before_table$after_table";
}

sub xref_total_summary {
  my ($self, $dbh, $title) = @_;

  my $sql = "
    SELECT
      db_name,
      COUNT(DISTINCT xref_id) AS xref_count,
      ensembl_object_type,
      COUNT(DISTINCT ensembl_id) AS ensembl_object_count,
      logic_name
    FROM
      analysis RIGHT OUTER JOIN
      object_xref USING (analysis_id) INNER JOIN
      xref USING (xref_id) INNER JOIN
      external_db USING (external_db_id)
    GROUP BY
      db_name,
      ensembl_object_type,
      logic_name
    ;
  ";

  my $sth = $dbh->prepare($sql);
  $sth->execute();

  my $columns = $sth->{NAME};
  my $results = $sth->fetchall_arrayref();

  return $self->format_table($title, $columns, $results);
}

sub xref_ontology_summary {
  my ($self, $dbh, $title) = @_;

  my $sql = "
    SELECT
      edb1.db_name,
      edb2.db_name AS source_db,
      linkage_type,
      COUNT(DISTINCT x1.xref_id) AS xref_count,
      ensembl_object_type,
      COUNT(DISTINCT ensembl_id) AS ensembl_object_count
    FROM
      xref x1 INNER JOIN
      object_xref ox ON x1.xref_id = ox.xref_id INNER JOIN
      ontology_xref ontx ON ox.object_xref_id = ontx.object_xref_id INNER JOIN
      xref x2 on ontx.source_xref_id = x2.xref_id INNER JOIN
      external_db edb1 on x1.external_db_id = edb1.external_db_id INNER JOIN
      external_db edb2 on x2.external_db_id = edb2.external_db_id
    GROUP BY
      edb1.db_name,
      edb2.db_name,
      linkage_type,
      ensembl_object_type
    ;
  ";

  my $sth = $dbh->prepare($sql);
  $sth->execute();

  my $columns = $sth->{NAME};
  my $results = $sth->fetchall_arrayref();

  return $self->format_table($title, $columns, $results);
}

1;
