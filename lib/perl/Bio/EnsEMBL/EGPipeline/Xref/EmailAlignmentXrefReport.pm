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

Bio::EnsEMBL::EGPipeline::Xref::EmailAlignmentXrefReport

=head1 DESCRIPTION

Run a few useful queries on the xrefs, for a given species.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::Xref::EmailAlignmentXrefReport;

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
  my $species     = $self->param_required('species');
  my $external_db = $self->param('external_db');
  my $logic_name  = $self->param('logic_name');

  my $dba = $self->get_DBAdaptor($self->param('db_type'));
  my $dbh = $dba->dbc->db_handle;

  my $reports = "The alignment xref pipeline for $species has completed.\n";
  $reports .= "Summaries are below; note that the last one includes pre-existing data.\n";

  $reports .= $self->xref_summary($dbh, $logic_name, "$external_db xrefs:");

  $reports .= $self->xref_total_summary($dbh, 'All xrefs, pre-existing and newly-added:');

  $self->param('text', $reports);
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

1;
