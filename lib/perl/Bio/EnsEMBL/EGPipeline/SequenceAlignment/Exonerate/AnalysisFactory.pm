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

package Bio::EnsEMBL::EGPipeline::SequenceAlignment::Exonerate::AnalysisFactory;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub param_defaults {
  return {};
}

sub run {
  my ($self) = @_;
  my $species = $self->param_required('species');
  my $analyses = $self->param_required('exonerate_analyses');
  my $analysis_name = $self->param_required('analysis_name');
  my $logic_name = $self->param_required('logic_name');
  my $db_backup_file = $self->param('db_backup_file');
  
  my $filtered_analyses = [];
  foreach my $analysis (@{$analyses}) {
    next unless $analysis_name eq $$analysis{'logic_name'};
    $$analysis{'logic_name'} = $logic_name;
    $$analysis{'species'} = $species;
    $$analysis{'db_backup_file'} = $db_backup_file if defined $db_backup_file;
    push @$filtered_analyses, $analysis;
  }
  $self->param('filtered_analyses', $filtered_analyses);
  
}

sub write_output {
  my ($self) = @_;
  
  $self->dataflow_output_id($self->param('filtered_analyses'), 2);
  
}

1;

