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

package Bio::EnsEMBL::EGPipeline::BlastAlignment::AnalysisUnmergeFactory;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub run {
  my ($self) = @_;
  my $program          = $self->param_required('program');
  my $external_db_name = $self->param_required('external_db_name');
  my $analysis_groups  = $self->param_required('analysis_groups');
  my $id_prefixes      = $self->param_required('id_prefixes');
  my $linked_tables    = $self->param_required('linked_tables');
  my $db_type          = $self->param_required('db_type');

  my $sql_commands = [];

  foreach my $table (@{$linked_tables}) {
    foreach my $group (keys %{$analysis_groups}) {
      foreach my $name (@{$$analysis_groups{$group}}) {
        if (! $self->has_data($db_type, $table, "$name\_$program")) {
          my $sql = $self->update_sql($table, "$group\_$program", "$name\_$program", $external_db_name, $$id_prefixes{$name});
          push @$sql_commands, {sql => $sql};
        }
      }
    }
  }

  $self->param('sql_commands', $sql_commands);
}

sub write_output {
  my ($self) = @_;

  $self->dataflow_output_id($self->param('sql_commands'), 2);
}

sub has_data {
  my ($self, $db_type, $table, $logic_name) = @_;

  my $dba = $self->get_DBAdaptor($db_type);
  my $dbh = $dba->dbc->db_handle();
  my $sql = "
    SELECT COUNT(*) FROM $table INNER JOIN analysis USING (analysis_id)
    WHERE logic_name = '$logic_name';
  ";
  my ($count) = $dbh->selectrow_array($sql);

  return $count ? 1 : 0;
}

sub update_sql {
  my ($self, $table, $from, $to, $external_db_name, $prefix) = @_;

  my $sql = "
    UPDATE
      $table        t INNER JOIN
      analysis     a1 USING (analysis_id) JOIN
      analysis     a2 JOIN
      external_db edb
    SET
      t.analysis_id    = a2.analysis_id,
      t.external_db_id = edb.external_db_id
    WHERE
      a1.logic_name = '$from' AND
      a2.logic_name = '$to' AND
      edb.db_name   = '$external_db_name' AND
      t.hit_name LIKE '$prefix\%';
  ";

  return $sql;
}

1;
