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

Bio::EnsEMBL::EGPipeline::PipeConfig::ProjectGenes_conf

=head1 DESCRIPTION

Project genes between two assemblies (expected to be similar,
e.g. different assembly versions or strains of the same species).

=cut

package Bio::EnsEMBL::EGPipeline::PipeConfig::ProjectGenes_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf');

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use Bio::EnsEMBL::Hive::Version 2.4;

use File::Spec::Functions qw(catdir);

sub default_options {
  my ($self) = @_;
  return {
    %{$self->SUPER::default_options},
    
    pipeline_name => 'project_genes',
    
    # Location of db with old assembly.
    from_db_url   => undef,
    
    # Location of db with new assembly.
    to_db_url     => undef,
    
    # The mapping process needs to write to the db with the old db. Since
    # that will usually be on a read-only archive server, by default make
    # a temporary copy for the pipeline to use.
    copy_from_db  => 1,
    
    # If analysis logic_names are not given, all genes will be projected.
    logic_name    => [],
    
    # Tidy and validate GFF3 result files.
    gt_exe        => 'gt',
    gff3_tidy     => $self->o('gt_exe').' gff3 -tidy -sort -retainids',
    gff3_validate => $self->o('gt_exe').' gff3validator',
  };
}

sub hive_meta_table {
  my ($self) = @_;

  return {
    %{$self->SUPER::hive_meta_table},
    'hive_use_param_stack'  => 1,
  };
}

sub pipeline_create_commands {
  my ($self) = @_;

  return [
    @{$self->SUPER::pipeline_create_commands},
    'mkdir -p '.$self->o('results_dir'),
  ];
}

sub pipeline_wide_parameters {
 my ($self) = @_;
 
 return {
   %{$self->SUPER::pipeline_wide_parameters},
   'copy_from_db' => $self->o('copy_from_db'),
 };
}

sub pipeline_analyses {
  my ($self) = @_;

  return [
    {
      -logic_name        => 'ProjectGenes',
      -module            => 'Bio::EnsEMBL::EGPipeline::ProjectGenes::ProjectGenes',
      -analysis_capacity => 5,
      -parameters        => {
                              from_db_url => $self->o('from_db_url'),
                              to_db_url   => $self->o('to_db_url'),
                              results_dir => $self->o('results_dir'),
                              logic_name  => $self->o('logic_name'),
                            },
      -input_ids         => [ {} ],
      -max_retry_count   => 0,
      -flow_into         => {
                              '2->A' => ['GFF3Tidy'],
                              'A->1' => ['EmailReport'],
                            },
      -rc_name           => '16Gb_mem',
    },

    {
      -logic_name        => 'GFF3Tidy',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -max_retry_count   => 0,
      -parameters        => {
                              cmd => $self->o('gff3_tidy').' #gff3_file# > #gff3_file#.tmp; '.
                                     'mv #gff3_file#.tmp #gff3_file#',
                            },
      -rc_name           => 'normal',
      -flow_into         => ['GFF3Validate'],
    },

    {
      -logic_name        => 'GFF3Validate',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -max_retry_count   => 0,
      -parameters        => {
                              cmd => $self->o('gff3_validate').' #gff3_file#',
                            },
      -rc_name           => 'normal',
    },

    {
      -logic_name        => 'EmailReport',
#      -module            => 'Bio::EnsEMBL::EGPipeline::ProjectGenes::EmailReport',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -max_retry_count   => 1,
      -parameters        => {
#                              email   => $self->o('email'),
#                              subject => 'Gene projection pipeline has completed for #species#',
                            },
      -rc_name           => 'normal',
    },

  ];
}

1;
