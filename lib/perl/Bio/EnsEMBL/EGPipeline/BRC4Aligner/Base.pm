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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::Base;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub store_align_cmds {
  my ($self, $cmds) = @_;
  
  my $query = "INSERT INTO align_cmds (sample_name, cmds, version) VALUES(?,?,?)";
  my $dbh = $self->hive_dbc->db_handle;
  my $sth = $dbh->prepare($query);

  my @values = map { $cmds->{$_} } qw(sample_name cmds version);
  $sth->execute(@values);
}

1;
