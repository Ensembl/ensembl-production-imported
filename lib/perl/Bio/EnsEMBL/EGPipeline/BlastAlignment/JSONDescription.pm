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

package Bio::EnsEMBL::EGPipeline::BlastAlignment::JSONDescription;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use File::Spec::Functions qw(catdir);

sub param_defaults {
  return {
    'db_type' => 'otherfeatures',
  };
}

sub run {
  my ($self) = @_;
  my $species     = $self->param_required('species');
  my $db_type     = $self->param_required('db_type');
  my $logic_name  = $self->param_required('logic_name')->[0];
  my $results_dir = $self->param_required('results_dir');

  my $dba = $self->get_DBAdaptor($db_type);
  my $aa = $dba->get_adaptor('Analysis');
  my $analysis = $aa->fetch_by_logic_name($logic_name);

  my $json_file = catdir($results_dir, "$species.$logic_name.json");

  open(my $fh, '>', $json_file) or $self->throw("Cannot open file $json_file: $!");

  print $fh
    '{
      "category":"protein_alignment",
      "version":"'.$analysis->db_version.'",
      "description":"'.$analysis->description.'",
      "source":"VectorBase"
    }';

  close($fh);

  $dba->dbc->disconnect_if_idle();
}

1;
