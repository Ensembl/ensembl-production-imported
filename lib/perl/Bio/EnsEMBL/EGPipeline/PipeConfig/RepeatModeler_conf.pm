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

package Bio::EnsEMBL::EGPipeline::PipeConfig::RepeatModeler_conf;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::Version 2.5;
use base ('Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf');

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

use File::Spec::Functions qw(catdir catfile rel2abs);

sub default_options {
  my ($self) = @_;
  return {
    %{$self->SUPER::default_options},

    pipeline_name => $self->tagged_pipeline_name('repeat_modeler'),
    work_dir => catdir($self->o('results_dir'), "work"),

    species => [],
    antispecies => [],
    division => [],
    run_all => 0,
    meta_filters => {},

    # Parameters for dumping and splitting Fasta DNA files
    max_seq_length          => 10_000_000,
    max_seq_length_per_file => $self->o('max_seq_length'),
    max_seqs_per_file       => 10_000,
    max_files_per_directory => 50,
    max_dirs_per_directory  => $self->o('max_files_per_directory'),
    min_slice_length        => 5000,

    # Program paths
    builddatabase_exe => 'BuildDatabase',
    repeatmodeler_exe => 'RepeatModeler',

    # Blast engine can be wublast or ncbi
    blast_engine => 'ncbi',

    # Cluster final library
    do_clustering       => 0,
    cdhit_est_exe       => $self->check_exe_in_cellar('cd-hit/4.6.8/bin/cd-hit-est'),
    
    # Filtering
    do_filtering => 0,
  };
}

# make use of init_pipeline.pl parameters
sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{$self->SUPER::pipeline_wide_parameters},
    
    # Splitting
    min_slice_length    =>  $self->o('min_slice_length'),
    
    # Dirs
    results_dir         => rel2abs($self->o('results_dir')),
    work_dir            => rel2abs($self->o('work_dir')),
    species_work_dir    => rel2abs(catdir($self->o('work_dir'), "#species#")),
    
    # Files
    repeat_lib          => catfile("#results_dir#", "#species#" . ".rm.lib"),
    repeat_lib_filtered => "#repeat_lib#" . ".filtered",
    transcripts_fasta   => catfile("#species_work_dir#", "#species#" . ".transcripts.fa"),
    translations_fasta  => catfile("#species_work_dir#", "#species#" . ".translations.fa"),
    repbase_species_lib  => catfile("#species_work_dir#", "#species#" . ".repbase.lib"),
    
    # Do special tasks
    do_filtering        => $self->o('do_filtering'),
    do_clustering       => $self->o('do_clustering'),
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
    'mkdir -p ' . $self->o('results_dir'),
    'mkdir -p ' . $self->o('work_dir'),
  ];
}

sub pipeline_analyses {
  my ($self) = @_;

  return [
    {
      -logic_name        => 'SpeciesFactory',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EGSpeciesFactory',
      -max_retry_count   => 1,
      -parameters        => {
                              species         => $self->o('species'),
                              antispecies     => $self->o('antispecies'),
                              division        => $self->o('division'),
                              run_all         => $self->o('run_all'),
                              meta_filters    => $self->o('meta_filters'),
                              chromosome_flow => 0,
                              regulation_flow => 0,
                              variation_flow  => 0,
                            },
      -input_ids         => [ {} ],
      -flow_into         => {
                              '2->A' => 'CreateLibrary',
                              'A->2' => 'Filtering',
                            },
      -meadow_type       => 'LOCAL',
    },

    {
      -logic_name        => 'CreateLibrary',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -max_retry_count   => 0,
      -parameters        => {
                              no_existing_lib => "#expr(not -e #repeat_lib#)expr#",
                            },
      -flow_into         => {
                              '1' => WHEN('#no_existing_lib#', 'CheckUsableLength'),
                            },
      -meadow_type       => 'LOCAL',
    },

    # The RepeatModeler needs a minimum length of sequences: don't try to run it below it
    {
      -logic_name        => 'CheckUsableLength',
      -module            => 'Bio::EnsEMBL::EGPipeline::DNAFeatures::CheckUsableLength',
      -parameters        => {
                              min_slice_length => $self->o('min_slice_length'),
                              repeat_modeler_length_threshold => 40_000,
                            },
      -flow_into         => {
                              '2' => 'DumpGenome',
                              '3' => 'NoModel'
                            },
      -analysis_capacity => 1,
      -max_retry_count   => 0,
      -batch_size        => 50,
      -meadow_type       => 'LOCAL',
    },

    # This is just to list all the species that do not have a model generated in the end
    {
      -logic_name        => 'NoModel',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters        =>
      {
        cmd => 'touch #repeat_lib#',
      },
      -analysis_capacity => 1,
      -max_retry_count   => 0,
      -batch_size        => 50,
      -meadow_type       => 'LOCAL',
    },

    # Same, when there is no filter needed
    {
      -logic_name        => 'NoFilter',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -analysis_capacity => 1,
      -max_retry_count   => 0,
      -batch_size        => 50,
      -meadow_type       => 'LOCAL',
    },

    {
      -logic_name        => 'DumpGenome',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::DumpGenome',
      -analysis_capacity => 5,
      -max_retry_count   => 1,
      -parameters        => {
                              genome_dir       => catdir("#work_dir#", '#species#'),
                              min_slice_length => $self->o('min_slice_length'),
                            },
      -rc_name           => 'default',
      -flow_into         => {
                              '1->A' => ['SplitDumpFiles'],
                              'A->1' => ['CheckConsensus'],
                            },
    },

    {
      -logic_name        => 'SplitDumpFiles',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::FastaSplit',
      -analysis_capacity => 5,
      -max_retry_count   => 0,
      -parameters        => {
                              fasta_file              => '#genome_file#',
                              max_seq_length_per_file => $self->o('max_seq_length_per_file'),
                              max_seqs_per_file       => $self->o('max_seqs_per_file'),
                              max_files_per_directory => $self->o('max_files_per_directory'),
                              max_dirs_per_directory  => $self->o('max_dirs_per_directory'),
                              unique_file_names       => 1,
                            },
      #-rc_name           => '8Gb_mem',
      -rc_name           => '16Gb_mem',
      -flow_into         => {'2' => ['BuildDatabase']},
    },

    {
      -logic_name        => 'BuildDatabase',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -analysis_capacity => 25,
      -max_retry_count   => 0,
      -batch_size        => 25,
      -parameters        =>
      {
        cmd => 'cd #work_dir#; '.
               'RM_DB=$(basename #split_file#); '.
               'mkdir $RM_DB; '.
               'cd $RM_DB; '.
               $self->o('builddatabase_exe').' -engine '.$self->o('blast_engine').' -name $RM_DB #split_file#',
      },
      -rc_name           => 'default',
      -flow_into         => ['RepeatModeler'],
    },

    {
      -logic_name        => 'RepeatModeler',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -analysis_capacity => 25,
      -max_retry_count   => 0,
      -parameters        =>
      {
        cmd => 'RM_DB=$(basename #split_file#); '.
               'cd #work_dir#/$RM_DB; '.
               $self->o('repeatmodeler_exe').' -pa 9 -engine '.$self->o('blast_engine').' -database $RM_DB',
      },
      -rc_name           => '16Gb_mem_8_cores',
      -flow_into         => {
                              '-1' => ['RepeatModeler_HighMem'],
                            },
    },

    {
      -logic_name        => 'RepeatModeler_HighMem',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -analysis_capacity => 25,
      -max_retry_count   => 0,
      -parameters        =>
      {
        cmd => 'RM_DB=$(basename #split_file#); '.
               'cd #work_dir#/$RM_DB; '.
               $self->o('repeatmodeler_exe').' -pa 9 -engine '.$self->o('blast_engine').' -database $RM_DB',
      },
      -rc_name           => '32Gb_mem_8_cores',
      -flow_into         => {
                              '-1' => ['RepeatModeler_HighMem_32'],
                            },
    },

    {
      -logic_name        => 'RepeatModeler_HighMem_32',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -analysis_capacity => 25,
      -max_retry_count   => 0,
      -parameters        =>
      {
        cmd => 'RM_DB=$(basename #split_file#); '.
               'cd #work_dir#/$RM_DB; '.
               $self->o('repeatmodeler_exe').' -pa 9 -engine '.$self->o('blast_engine').' -database $RM_DB',
      },
      -rc_name           => '64Gb_mem_8_cores',
    },

    {
      -logic_name        => 'CheckConsensus',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -max_retry_count   => 0,
      -parameters        =>
      {
        consensus_file => '#work_dir#/#species#*/RM_*/consensi.fa.classified',
        has_consensus => '#expr(glob(#consensus_file#))expr#',
      },
      -flow_into         => {
                              '1' => WHEN("#has_consensus#", 'MergeResults', ELSE({ 'NoModel' => { species => '#species#', reason => 'no_consensus' } })),
                            },
      -analysis_capacity => 1,
      -meadow_type       => 'LOCAL',
    },

    # MERGE AND CLUSTER
    {
      -logic_name        => 'MergeResults',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -analysis_capacity => 5,
      -max_retry_count   => 0,
      -parameters        =>
      {
        out_file => "#repeat_lib#",
        cmd => 'cat #work_dir#/#species#*/RM_*/consensi.fa.classified > #out_file#',
      },
      -rc_name           => 'normal',
      -flow_into         => WHEN('#do_clustering#' => 'ClusterLib'),
    },
    
    {
      -logic_name        => 'ClusterLib',
      -module            => 'Bio::EnsEMBL::EGPipeline::DNAFeatures::ClusterRepeatLib',
      -max_retry_count   => 0,
      -parameters        =>
      {
        lib_file      => '#repeat_lib#',
        out_file      => '#repeat_lib#',
        wd_path       => '#species_work_dir#/repeat_clustering',

        cdhit_est_exe            => $self->o('cdhit_est_exe'),
        cdhit_identity_threshold => 0.8,
        cdhit_alignment_coverage => 0.9,
        cdhit_word_len           => 5,
        cdhit_num_threads        => 4,
      },
      -rc_name           => 'normal',
    },

    # FILTER
    {
      -logic_name        => 'Filtering',
      -module            => 'Bio::EnsEMBL::EGPipeline::DNAFeatures::HasGenes',
      -parameters        =>
      {
        need_filtering      => '#expr(-s #repeat_lib# and not -e #repeat_lib_filtered#)expr#',
      },
      -flow_into         => {
                              '2' => WHEN('#do_filtering# and #need_filtering#', 'Dumps', ELSE({ 'NoFilter' => { species => '#species#' } })),
                            },
      -meadow_type       => 'LOCAL',
    },

    {
      -logic_name        => 'Dumps',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters        =>
      {
        cmd => 'mkdir -p #species_work_dir#',
      },
      -flow_into         => {
        '1->A' => [
          'DumpTranscripts',
          'DumpTranslations',
          'DumpSpeciesRepbase',
        ],
        'A->1' => 'RepbaseTranslationsToIgnore',
      },
      -meadow_type       => 'LOCAL',
    },
    {
      -logic_name        => 'DumpTranscripts',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::DumpTranscriptome',
      -analysis_capacity => 5,
      -parameters        =>
      {
        transcriptome_file      => '#transcripts_fasta#',
        overwrite       => 1,
      },
      -rc_name           => 'normal',
    },
    {
      -logic_name        => 'DumpTranslations',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::DumpProteome',
      -analysis_capacity => 5,
      -parameters        =>
      {
        proteome_file   => '#translations_fasta#',
        overwrite       => 1,
      },
      -rc_name           => 'normal',
    },
    {
      -logic_name        => 'DumpSpeciesRepbase',
      -module            => 'Bio::EnsEMBL::EGPipeline::DNAFeatures::DumpSpeciesRepbase',
      -analysis_capacity => 5,
      -parameters        =>
      {
        species_repeat_lib => '#repbase_species_lib#',
      },
      -max_retry_count   => 0,
      -rc_name           => 'normal',
    },

    {
      -logic_name        => 'RepbaseTranslationsToIgnore',
      -module            => 'Bio::EnsEMBL::EGPipeline::DNAFeatures::RepbaseTranslationsToIgnore',
      -analysis_capacity => 5,
      -parameters        =>
      {
        transcripts_file   => '#transcripts_fasta#',
        translations_file  => '#translations_fasta#',
        species_repeat_lib => '#repbase_species_lib#',
        threads            => 1,
      },
      -max_retry_count   => 0,
      -flow_into         => {
        1 => 'FilterCustomLib',
      },
      -rc_name           => 'normal',
    },
    {
      -logic_name        => 'FilterCustomLib',
      -module            => 'Bio::EnsEMBL::EGPipeline::DNAFeatures::FilterCustomLib',
      -analysis_capacity => 5,
      -parameters        =>
      {
        transcripts_file   => '#transcripts_fasta#',
        lib_to_filter      => "#repeat_lib#",
        filtered_lib       => "#repeat_lib_filtered#",
        threads            => 1,
      },
      -max_retry_count   => 0,
      -rc_name           => 'normal',
    },
  ];
}

sub resource_classes {
  my ($self) = @_;
  
  return {
    %{$self->SUPER::resource_classes},
    '8Gb_mem_8_cores'  => {'LSF' => '-q ' . $self->o('queue_name') . ' -M 8000 -n 8 -R "rusage[mem=8000,tmp=4000]"'},
    '16Gb_mem_8_cores' => {'LSF' => '-q ' . $self->o('queue_name') . ' -M 16000 -n 8 -R "rusage[mem=16000,tmp=4000]"'},
    '32Gb_mem_8_cores' => {'LSF' => '-q ' . $self->o('queue_name') . ' -M 32000 -n 8 -R "rusage[mem=32000,tmp=4000]"'},
    '64Gb_mem_8_cores' => {'LSF' => '-q ' . $self->o('queue_name') . '  -M 64000 -n 8 -R "rusage[mem=64000,tmp=4000]"'},
  }
}

1;
