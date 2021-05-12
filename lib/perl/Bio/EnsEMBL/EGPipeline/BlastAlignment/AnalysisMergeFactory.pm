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

package Bio::EnsEMBL::EGPipeline::BlastAlignment::AnalysisMergeFactory;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub run {
  my ($self) = @_;
  my $program          = $self->param_required('program');
  my $external_db_name = $self->param_required('external_db_name');
  my $analysis_groups  = $self->param_required('analysis_groups');
  my $logic_name       = $self->param_required('logic_name');
  my $linked_tables    = $self->param_required('linked_tables');
  
  (my $group = $logic_name) =~ s/_$program$//;
  
  my $sql_commands = [];
  
  foreach my $table (@{$linked_tables}) {
    foreach my $name (@{$$analysis_groups{$group}}) {
      my $sql = $self->update_sql($table, "$name\_$program", $logic_name, $external_db_name);
      push @$sql_commands, {sql => $sql};
    }
  }
  
  $self->param('sql_commands', $sql_commands);
}

sub write_output {
  my ($self) = @_;
  
  $self->dataflow_output_id($self->param('sql_commands'), 2);
}

sub update_sql {
  my ($self, $table, $from, $to, $external_db_name) = @_;
  
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
      edb.db_name   = '$external_db_name';
  ";
  
  return $sql;
}

1;

