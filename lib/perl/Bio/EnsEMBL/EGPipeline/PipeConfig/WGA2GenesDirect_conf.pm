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

Bio::EnsEMBL::EGPipeline::PipeConfig::WGA2GenesDirect_conf

=head1 DESCRIPTION

Project transripts using compara produced lastz mappings.

=cut

package Bio::EnsEMBL::EGPipeline::PipeConfig::WGA2GenesDirect_conf;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::Version 2.5;
use base ('Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf');

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub default_options {
  my ($self) = @_;
  return {
    %{$self->SUPER::default_options},

    pipeline_name => $self->tagged_pipeline_name('wga2genes', $self->o('ensembl_release')),

    # run options 
    'projection_exonerate_padding' => 50,
    'method_link_type'             => 'LASTZ_NET',
    'result_force_rewrite'         => 0,
    'result_clone_mode'            => 'clone', # or dna_db

    # db options
    'driver'  => 'mysql',
    'compara_pass' => undef,    
    'source_pass' => undef,    
    'target_pass' => undef,    

    # compara db with lastz 
    'compara_db' => {
      -host   => $self->o('compara_host'),
      -port   => $self->o('compara_port'),
      -user   => $self->o('compara_user'),
      -pass   => $self->o('compara_pass'),
      -dbname => $self->o('compara_dbname'),
    },

    # source species (using same dna and transcript db)
    'source_db' => {
      -host   => $self->o('source_host'),
      -port   => $self->o('source_port'),
      -user   => $self->o('source_user'),
      -pass   => $self->o('source_pass'),
      -dbname => $self->o('source_dbname'),
    },

    # target species
    'target_db' => {
      -host   => $self->o('target_host'),
      -port   => $self->o('target_port'),
      -user   => $self->o('target_user'),
      -pass   => $self->o('target_pass'),
      -dbname => $self->o('target_dbname'),
    },

    # results / projection
    'result_db' => {
      -host   => $self->o('result_host'),
      -port   => $self->o('result_port'),
      -user   => $self->o('result_user'),
      -pass   => $self->o('result_pass'),
      -dbname => $self->o('result_dbname'),
    },
  }; # end return
} # end default_options

sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{$self->SUPER::pipeline_wide_parameters},
  };
}

# Force an automatic loading of the registry in all workers.
sub beekeeper_extra_cmdline_options {
  my ($self) = @_;

  my $options = join(' ',
    $self->SUPER::beekeeper_extra_cmdline_options,
    "-reg_conf ".$self->o('registry')
  );
  
  return $options;
}

# Ensures that species output parameter gets propagated implicitly.
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
      #'mkdir -p '.$self->o('results_dir'),
    ];
} # end pipeline_create_commands


sub pipeline_analyses {
  my ($self) = @_;

  return [
    {
      -logic_name => 'create_output_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
          source_db => $self->o('target_db'),
          target_db => $self->o('result_db'),
          create_type => $self->o('result_clone_mode'),
          _lock_tables => 'false',
          force_drop => $self->o('result_force_rewrite'),
        },
      -input_ids => [{}],
      -flow_into => {
          1 => ['generate_projection_ids'],
      },
    },

    {
      -logic_name => 'generate_projection_ids',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
          target_db => $self->o('source_db'),
          iid_type => 'feature_id',
          feature_type => 'transcript',
          feature_restriction => 'projection',
          biotypes => {
            'protein_coding' => 1,
          },
          batch_size => 200,
      },
      -flow_into => {
          2 => ['project_transcripts'],
      },
      -rc_name    => 'default',
    },

    {
      -logic_name => 'project_transcripts',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveWGA2GenesDirect',
      -parameters => {
          logic_name => 'project_transcripts',
          module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveWGA2GenesDirect',
          source_transcript_db => $self->o('source_db'),
          source_dna_db => $self->o('source_db'),
          target_transcript_db => $self->o('result_db'),
          target_dna_db => $self->o('target_db'),
          compara_db => $self->o('compara_db'),
          method_link_type => $self->o('method_link_type'),
          MAX_EXON_READTHROUGH_DIST => 15,
          TRANSCRIPT_FILTER => {
            OBJECT     => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
            PARAMETERS => {
              -coverage => 50,
              -percent_id => 50,
            },
          },
          iid_type => 'feature_id',
          feature_type => 'transcript',
          calculate_coverage_and_pid => 1,
          max_internal_stops => 1,
          timer => '10m', # value in the form of 1h20m or 2H05M
      },
      -flow_into => {
          -3 => ['failed_projection_jobs'],
       },
      -rc_name       => 'default',
      -hive_capacity => 20,
    },

    {
      -logic_name => 'failed_projection_jobs',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -parameters => {
      },
      -rc_name          => 'default',
      -can_be_empty  => 1,
    },
  ];
} # end pipeline analyses

1;
