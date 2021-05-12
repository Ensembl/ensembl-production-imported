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

package Bio::EnsEMBL::EGPipeline::BlastAlignment::AnalysisSetupFactory;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub run {
  my ($self) = @_;
  my $program         = $self->param_required('program');
  my $analysis_groups = $self->param_required('analysis_groups');
  my $analyses        = $self->param_required('analyses');
  
  my $merged_analyses = [];
  
  foreach my $analysis_group (keys %{$analysis_groups}) {
    foreach my $analysis (@{$analyses}) {
      my %merged_analysis = %$analysis;
      $merged_analysis{'logic_name'} = "$analysis_group\_$program";
      push @$merged_analyses, \%merged_analysis;
    }
  }
  
  $self->param('merged_analyses', $merged_analyses);
}

sub write_output {
  my ($self) = @_;
  
  $self->dataflow_output_id($self->param('merged_analyses'), 2);
}

1;
