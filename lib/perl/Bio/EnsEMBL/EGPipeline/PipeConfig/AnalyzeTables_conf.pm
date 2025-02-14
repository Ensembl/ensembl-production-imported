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

Bio::EnsEMBL::EGPipeline::PipeConfig::AnalyzeTables_conf

=head1 DESCRIPTION

Analyze (or optionally optimize) all tables in a database.

=cut

package Bio::EnsEMBL::EGPipeline::PipeConfig::AnalyzeTables_conf;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::Version v2.3;
use base ('Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf');

sub default_options {
  my ($self) = @_;
  return {
    %{$self->SUPER::default_options},
    
    pipeline_name => 'analyze_tables_'.$self->o('ensembl_release'),
    
    species      => [],
    division     => [],
    run_all      => 0,
    antispecies  => [],
    meta_filters => {},
    
    optimize_tables => 0,
  };
}

sub beekeeper_extra_cmdline_options {
  my $self = shift;
  return "-reg_conf ".$self->o("registry");
}

sub pipeline_analyses {
  my ($self) = @_;
  
  return [
    {
      -logic_name      => 'SpeciesFactory',
      -module          => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EGSpeciesFactory',
      -parameters      => {
                            species         => $self->o('species'),
                            division        => $self->o('division'),
                            run_all         => $self->o('run_all'),
                            antispecies     => $self->o('antispecies'),
                            meta_filters    => $self->o('meta_filters'),
                            core_flow       => 2,
                            chromosome_flow => 0,
                            regulation_flow => 3,
                            variation_flow  => 4,
                          },
      -input_ids       => [ {} ],
      -max_retry_count => 1,
      -flow_into       => {
                            '2' => ['AnalyzeTablesCore', 'AnalyzeTablesOtherFeatures'],
                            '3' => ['AnalyzeTablesRegulation'],
                            '4' => ['AnalyzeTablesVariation'],
                          },
      -meadow_type     => 'LOCAL',
    },

    {
      -logic_name        => 'AnalyzeTablesCore',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::AnalyzeTables',
      -parameters        => {
                              optimize_tables => $self->o('optimize_tables'),
                              db_types        => ['core'],
                            },
      -max_retry_count   => 1,
      -analysis_capacity => 10,
      -rc_name           => 'normal',
    },
    
    {
      -logic_name        => 'AnalyzeTablesOtherFeatures',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::AnalyzeTables',
      -parameters        => {
                              optimize_tables => $self->o('optimize_tables'),
                              db_types        => ['otherfeatures'],
                            },
      -max_retry_count   => 1,
      -analysis_capacity => 10,
      -rc_name           => 'normal',
    },
    
    {
      -logic_name        => 'AnalyzeTablesRegulation',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::AnalyzeTables',
      -parameters        => {
                              optimize_tables => $self->o('optimize_tables'),
                              db_types        => ['funcgen'],
                            },
      -max_retry_count   => 1,
      -analysis_capacity => 10,
      -rc_name           => 'normal',
    },
    
    {
      -logic_name        => 'AnalyzeTablesVariation',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::AnalyzeTables',
      -parameters        => {
                              optimize_tables => $self->o('optimize_tables'),
                              db_types        => ['variation'],
                            },
      -max_retry_count   => 1,
      -analysis_capacity => 10,
      -rc_name           => 'normal',
    },
  ];
}

1;
