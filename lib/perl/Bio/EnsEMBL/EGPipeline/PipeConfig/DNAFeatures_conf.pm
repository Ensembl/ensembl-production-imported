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
Bio::EnsEMBL::EGPipeline::PipeConfig::DNAFeatures_conf

=head1 DESCRIPTION

Configuration for running the DNA Features pipeline, which
primarily adds repeat features to a core database.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::PipeConfig::DNAFeatures_conf;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use Bio::EnsEMBL::Hive::Version 2.4;

use base ('Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf');
use File::Spec::Functions qw(catdir catfile rel2abs);

sub default_options {
  my ($self) = @_;

  return {
    %{$self->SUPER::default_options},

    parameters_json  => undef,

    pipeline_name => $self->tagged_pipeline_name('dna_features', $self->o('ensembl_release')),
    report_dir    => catdir($self->o('pipeline_dir'), 'reports'),

    species      => [],
    antispecies  => [],
    taxons       => [],
    antitaxons   => [],
    division     => [],
    run_all      => 0,
    meta_filters => {},
    
    # Production can be given by the registry, or by this
    production_db => undef,

    # Default parameters for dumping and splitting Fasta DNA files.
    max_seq_length          => 1000000,
    max_seq_length_per_file => $self->o('max_seq_length'),
    max_seqs_per_file       => undef,
    max_files_per_directory => 50,
    max_dirs_per_directory  => $self->o('max_files_per_directory'),
    
    # Dust and TRF can handle large files; this size should mean that
    # jobs take a few minutes.
    dust_trf_max_seq_length    => 100000000,
    dust_trf_max_seqs_per_file =>     10000,
    
    # Values >100 are not recommended, because you tend to overload
    # the mysql server with connections.
    max_hive_capacity => 100,

    dust_exe         => 'dustmasker',
    trf_exe          => 'trf',
    repeatmasker_exe => 'RepeatMasker',

    dust              => 1,
    trf               => 1,
    repeatmasker      => 1,
    redatrepeatmasker => 0,

    # By default, run RepeatMasker with repbase library and exclude
    # low-complexity annotations. By explicitly turning on the GC
    # calculations, the results are made consistent regardless of the
    # number of sequences in the input file. For repbase, a species
    # parameter is added when the program is called within the
    # pipeline. The sensitivity of the search, including which engine
    # is used, is also added within the pipeline.
    always_use_repbase         => 0,
    repeatmasker_library       => {},
    redat_repeatmasker_library => $self->private_conf('REDAT_REPEATMASKER_LIBRARY_PATH'),
    repeatmasker_sensitivity   => {},
    repeatmasker_logic_name    => {},
    repeatmasker_parameters    => ' -nolow -gccalc ',
    repeatmasker_cache         => catdir($self->o('pipeline_dir'), 'cache'),
    repeatmasker_timer         => '18H',
    repeatmasker_resource_class => 'normal',
    
    # Override species name to use with repbase
    repeatmasker_repbase_species => '',
    
    # Instead, use the species classification to find the closest species in Repbase
    guess_repbase_species => 0,

    # The ensembl-analysis Dust and TRF modules take a parameters hash
    # which is parsed, rather than requiring explicit command line
    # options.  It's generally not necessary to override default
    # values, but below are examples of the syntax for dust and trf,
    # showing the current defaults.  (See the help for those programs
    # for parameter descriptions.)

    # dust_parameters_hash => {
    #   'MASKING_THRESHOLD' => 20,
    #   'WORD_SIZE'         => 3,
    #   'WINDOW_SIZE'       => 64,
    #   'SPLIT_LENGTH'      => 50000,
    # },

    # trf_parameters_hash => {
    #   'MATCH'      => 2,
    #   'MISMATCH'   => 5,
    #   'DELTA'      => 7,
    #   'PM'         => 80,
    #   'PI'         => 10,
    #   'MINSCORE'   => 40,
    #   'MAX_PERIOD' => 500,
    # },

    dust_parameters_hash => {},
    trf_parameters_hash  => {},

    dna_analyses =>
    [
      {
        'logic_name'      => 'dust',
        'program'         => 'dustmasker',
        'program_file'    => $self->o('dust_exe'),
        'module'          => 'Bio::EnsEMBL::Analysis::Runnable::DustMasker',
        'gff_source'      => 'dust',
        'gff_feature'     => 'low_complexity_region',
        'linked_tables'   => ['repeat_feature'],
      },
      {
        'logic_name'      => 'trf',
        'program'         => 'trf',
        'program_version' => '4.0',
        'program_file'    => $self->o('trf_exe'),
        'module'          => 'Bio::EnsEMBL::Analysis::Runnable::TRF',
        'gff_source'      => 'trf',
        'gff_feature'     => 'tandem_repeat',
        'linked_tables'   => ['repeat_feature'],
      },
      {
        'logic_name'      => 'repeatmask_redat',
        'db'              => 'redat',
        'program'         => 'RepeatMasker',
        'program_version' => '4.0.5',
        'program_file'    => $self->o('repeatmasker_exe'),
        'parameters'      => $self->o('repeatmasker_parameters'),
        'module'          => 'Bio::EnsEMBL::Analysis::Runnable::REdatRepeatMasker',
        'gff_source'      => 'repeatmasker',
        'gff_feature'     => 'repeat_region',
        'linked_tables'   => ['repeat_feature'],
        'timer'           => '16H',
      },
      {
        'logic_name'      => 'repeatmask_repbase',
        'db'              => 'repbase',
        'program'         => 'RepeatMasker',
        'program_version' => '4.0.5',
        'program_file'    => $self->o('repeatmasker_exe'),
        'parameters'      => $self->o('repeatmasker_parameters'),
        'module'          => 'Bio::EnsEMBL::Analysis::Runnable::RepeatMasker',
        'gff_source'      => 'repeatmasker',
        'gff_feature'     => 'repeat_region',
        'linked_tables'   => ['repeat_feature'],
        'timer'           => '16H',
        'repbase_species'     => $self->o('repeatmasker_repbase_species')
      },
      {
        'logic_name'      => 'repeatmask_customlib',
        'db'              => 'custom',
        'program'         => 'RepeatMasker',
        'program_version' => '4.0.5',
        'program_file'    => $self->o('repeatmasker_exe'),
        'parameters'      => $self->o('repeatmasker_parameters'),
        'module'          => 'Bio::EnsEMBL::Analysis::Runnable::RepeatMasker',
        'gff_source'      => 'repeatmasker',
        'gff_feature'     => 'repeat_region',
        'linked_tables'   => ['repeat_feature'],
      },
    ],

    # Remove existing DNA features; if => 0 then existing analyses and
    # their features will remain, with the logic_name suffixed by
    # '_bkp'.
    delete_existing => 1,

    # Retrieve analysis descriptions from the production database; the
    # supplied registry file will need the relevant server details.
    production_lookup => 1,

    # By default, an email is sent for each species when the pipeline
    # is complete, showing the breakdown of repeat coverage.
    email_report => 1,
  };
}

# Force an automatic loading of the registry in all workers. WHY
sub beekeeper_extra_cmdline_options {
  my ($self) = @_;

  my $options = join(' ',
    $self->SUPER::beekeeper_extra_cmdline_options,
    "-reg_conf ".$self->o('registry')
  );

  return $options;
}

# Ensures that species output parameter gets propagated implicitly. WHY
sub hive_meta_table {
  my ($self) = @_;

  return {
    %{$self->SUPER::hive_meta_table},
    'hive_use_param_stack'  => 1,
  };
}


# Ensure the repeatmasker_cache directory exists
sub pipeline_create_commands {
  my ($self) = @_;

  return [
    @{$self->SUPER::pipeline_create_commands},
    'mkdir -p '.$self->o('repeatmasker_cache'),
  ];
}

sub pipeline_wide_parameters {
 my ($self) = @_;

 return {
   %{$self->SUPER::pipeline_wide_parameters},
    parameters_json  => $self->o('parameters_json'),
    work_dir                        => rel2abs($self->o('pipeline_dir')),
   'dust'                           => $self->o('dust'),
   'trf'                            => $self->o('trf'),
   'repeatmasker'                   => $self->o('repeatmasker'),
   'redatrepeatmasker'              => $self->o('redatrepeatmasker'),
   'email_report'                   => $self->o('email_report'),
   'repeatmasker_timer'             => $self->o('repeatmasker_timer'),
   'repeatmasker_repbase_species'   => $self->o('repeatmasker_repbase_species'),
   'repeatmasker_library'           => $self->o('repeatmasker_library'),
   'repeatmasker_sensitivity'       => $self->o('repeatmasker_sensitivity'),
   'repeatmasker_logic_name'        => $self->o('repeatmasker_logic_name'),
   'repeatmasker_resource_class'    => $self->o('repeatmasker_resource_class'),
 };
}

sub pipeline_analyses {
  my ($self) = @_;

  return [
    {
      -logic_name        => 'ParametersFromJson',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::ParametersFromJson',
      -analysis_capacity => 1,
      -input_ids         => [ {} ],
      -max_retry_count   => 0,
      -parameters        => {
                              parameters_json  => $self->o('parameters_json'),
                            },
      -flow_into         => {
                              '1' => ['DBFactory'],
                            },
      -meadow_type       => 'LOCAL',
    },
    {
      -logic_name        => 'DBFactory',
      -module            => 'Bio::EnsEMBL::Production::Pipeline::Common::DbFactory',
      -max_retry_count   => 1,
      -parameters        => {
                              species         => $self->o('species'),
                              antispecies     => $self->o('antispecies'),
                              taxons          => $self->o('taxons'),
                              antitaxons      => $self->o('antitaxons'),
                              division        => $self->o('division'),
                              run_all         => $self->o('run_all'),
                              meta_filters    => $self->o('meta_filters'),
                            },
      -flow_into         => {
                              '2->A' => ['BackupDatabase'],
                              'A->2' => ['UpdateMetadata'],
                            },
      -meadow_type       => 'LOCAL',
    },

    {
      -logic_name        => 'BackupDatabase',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::DatabaseDumper',
      -analysis_capacity => 5,
      -max_retry_count   => 1,
      -parameters        => {
                              output_file => catdir('#work_dir#', '#dbname#_bkp.sql.gz'),
                              table_list  => ['analysis', 'analysis_description', 'repeat_feature', 'repeat_consensus'],
                            },
      -rc_name           => 'normal',
      -flow_into         => {
                              '1->A' => ['DNAAnalysisFactory'],
                              'A->1' => ['DbAwareSpeciesFactory'],
                            },
    },

    {
      -logic_name        => 'DNAAnalysisFactory',
      -module            => 'Bio::EnsEMBL::EGPipeline::DNAFeatures::DNAAnalysisFactory',
      -max_retry_count   => 0,
      -analysis_capacity => 10,
      -batch_size        => 10,
      -parameters        => {
                              dust               => $self->o('dust'),
                              trf                => $self->o('trf'),
                              repeatmasker       => $self->o('repeatmasker'),
                              redatrepeatmasker  => $self->o('redatrepeatmasker'),

                              dna_analyses       => $self->o('dna_analyses'),
                              max_seq_length     => $self->o('max_seq_length'),
                              always_use_repbase => $self->o('always_use_repbase'),
                              rd_rm_library      => $self->o('redat_repeatmasker_library'),
                              rm_library         => '#repeatmasker_library#',
                              rm_sensitivity     => '#repeatmasker_sensitivity#',
                              rm_logic_name      => '#repeatmasker_logic_name#',
                              pipeline_dir       => '#work_dir#',
                              db_backup_file     => catdir('#work_dir#', '#dbname#_bkp.sql.gz'),
                              guess_repbase_species => $self->o('guess_repbase_species'),
                            },
      -rc_name           => 'normal',
      -flow_into         => {
                              '2->A' => ['AnalysisSetup'],
                              'A->1' => ['DeleteRepeatConsensus'],
                            },
    },

    {
      -logic_name        => 'AnalysisSetup',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::AnalysisSetup',
      -max_retry_count   => 0,
      -analysis_capacity => 10,
      -batch_size        => 10,
      -parameters        => {
                              db_backup_required => 1,
                              delete_existing    => $self->o('delete_existing'),
                              production_lookup  => $self->o('production_lookup'),
                              production_db      => $self->o('production_db'),
                            },
      -meadow_type       => 'LOCAL',
    },

    {
      -logic_name        => 'DeleteRepeatConsensus',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::SqlCmd',
      -max_retry_count   => 1,
      -analysis_capacity => 10,
      -batch_size        => 10,
      -parameters        => {
                             description => 'Clean up repeat_consensus features',
                             sql => [
                                'DELETE rc.* FROM '.
                                'repeat_consensus rc LEFT OUTER JOIN '.
                                'repeat_feature rf USING (repeat_consensus_id) '.
                                'WHERE rf.repeat_consensus_id IS NULL'
                              ]
                            },
      -meadow_type       => 'LOCAL',
    },

    {
      -logic_name        => 'DbAwareSpeciesFactory',
      -module            => 'Bio::EnsEMBL::Production::Pipeline::Common::DbAwareSpeciesFactory',
      -max_retry_count   => 1,
      -analysis_capacity => 2,
      -parameters        => {},
      -rc_name           => 'normal',
      -flow_into         => {
                              '2' => ['DumpGenome'],
                            },
    },

    {
      -logic_name        => 'DumpGenome',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::DumpGenome',
      -analysis_capacity => 5,
      -parameters        => {
                              genome_dir => catdir('#work_dir#', '#species#'),
                            },
      -rc_name           => 'normal',
      -flow_into         => [
                              WHEN('#dust# || #trf#' => ['SplitDumpFiles_1']),
                              WHEN('#repeatmasker# || #redatrepeatmasker#'  => ['RepeatMaskerWD']),
                            ],
    },

    {
      -logic_name        => 'SplitDumpFiles_1',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::FastaSplit',
      -analysis_capacity => 25,
      -parameters        => {
                              fasta_file              => '#genome_file#',
                              max_seq_length_per_file => $self->o('dust_trf_max_seq_length'),
                              max_seqs_per_file       => $self->o('dust_trf_max_seqs_per_file'),
                              max_files_per_directory => $self->o('max_files_per_directory'),
                              max_dirs_per_directory  => $self->o('max_dirs_per_directory'),
                              out_dir                 => catdir('#work_dir#', '#species#', 'dust_trf'),
                              file_varname            => 'queryfile',
                            },
      -rc_name           => '8Gb_mem',
      -flow_into         => {
                              '2' => [
                                WHEN('#dust#' => ['Dust']),
                                WHEN('#trf#'  => ['TRF']),
                              ],
                            },
    },

    {
      -logic_name        => 'RepeatMaskerWD',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters         => {
                               repeat_masker_wd => catdir('#work_dir#', '#species#', 'repeat_masker_wd'),
                               cmd => 'mkdir -p #repeat_masker_wd#',
      },
      -rc_name           => 'normal',
      -max_retry_count => 0,
      -flow_into         => {
                              '1->A' => ['SplitDumpFiles_2'],
                              'A->1' => ['RepeatMaskerWDCleanUp'],
                            },
    },

    {
      -logic_name        => 'RepeatMaskerWDCleanUp',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters         => {
                               repeat_masker_wd => catdir('#work_dir#', '#species#', 'repeat_masker_wd'),
                               cmd => 'find #repeat_masker_wd# -type f -print0 | xargs -r -0 rm',
      },
      -rc_name           => 'normal',
      -max_retry_count => 0,
    },

    {
      -logic_name        => 'SplitDumpFiles_2',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::FastaSplit',
      -analysis_capacity => 25,
      -parameters        => {
                              fasta_file              => '#genome_file#',
                              max_seq_length_per_file => $self->o('max_seq_length_per_file'),
                              max_seqs_per_file       => $self->o('max_seqs_per_file'),
                              max_files_per_directory => $self->o('max_files_per_directory'),
                              max_dirs_per_directory  => $self->o('max_dirs_per_directory'),
                              out_dir                 => catdir('#work_dir#', '#species#', 'repeatmasker'),
                              file_varname            => 'queryfile',
                            },
      -rc_name           => '8Gb_mem',
      -flow_into         => {
                              '2' => [
                                WHEN('#repeatmasker#'      => ['RepeatMaskerFactory']),
                                WHEN('#redatrepeatmasker#' => ['REdatRepeatMaskerFactory']),
                              ],
                            },
    },

    {
      -logic_name        => 'Dust',
      -module            => 'Bio::EnsEMBL::EGPipeline::DNAFeatures::DustMasker',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -max_retry_count   => 1,
      -batch_size        => 100,
      -parameters        => {
                              logic_name      => 'dust',
                              parameters_hash => $self->o('dust_parameters_hash'),
                            },
      -rc_name           => '16Gb_mem',
    },

    {
      -logic_name        => 'TRF',
      -module            => 'Bio::EnsEMBL::EGPipeline::DNAFeatures::TRF',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -max_retry_count   => 1,
      -batch_size        => 10,
      -parameters        => {
                              logic_name      => 'trf',
                              parameters_hash => $self->o('trf_parameters_hash'),
                            },
      -rc_name           => 'normal',
    },

    {
      -logic_name        => 'RepeatMaskerFactory',
      -module            => 'Bio::EnsEMBL::EGPipeline::DNAFeatures::RepeatMaskerFactory',
      -max_retry_count   => 1,
      -analysis_capacity => 10,
      -batch_size        => 10,
      -parameters        => {
                              always_use_repbase => $self->o('always_use_repbase'),
                              rm_library         => '#repeatmasker_library#',
                              rm_logic_name      => '#repeatmasker_logic_name#',
                              max_seq_length     => $self->o('max_seq_length'),
                            },
      -rc_name           => '16Gb_mem',
      -flow_into         => ['RepeatMasker'],
    },

    {
      -logic_name        => 'RepeatMasker',
      -module            => 'Bio::EnsEMBL::EGPipeline::DNAFeatures::RepeatMasker',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -max_retry_count   => 1,
      -parameters        => {
                              repeatmasker_cache => $self->o('repeatmasker_cache'),
                              timer              => $self->o('repeatmasker_timer'),
                              workdir            => catdir('#work_dir#', '#species#', 'repeat_masker_wd'),
                            },
      -rc_name           => $self->o('repeatmasker_resource_class'),
    },

    {
      -logic_name        => 'REdatRepeatMaskerFactory',
      -module            => 'Bio::EnsEMBL::EGPipeline::DNAFeatures::REdatRepeatMaskerFactory',
      -max_retry_count   => 1,
      -analysis_capacity => 10,
      -batch_size        => 10,
      -parameters        => {
                              rd_rm_library        => $self->o('redat_repeatmasker_library'),
                              logic_name           => 'repeatmask_redat',
                              max_seq_length       => $self->o('max_seq_length'),
                            },
      -rc_name           => '8Gb_mem',
      -flow_into         => ['REdatRepeatMasker'],
    },

    {
      -logic_name        => 'REdatRepeatMasker',
      -module            => 'Bio::EnsEMBL::EGPipeline::DNAFeatures::REdatRepeatMasker',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -max_retry_count   => 1,
      -parameters        => {
                              repeatmasker_cache => $self->o('repeatmasker_cache'),
                              timer              => $self->o('repeatmasker_timer'),
                              workdir            => catdir('#work_dir#', '#species#', 'repeat_masker_wd'),
                            },
      -rc_name           => 'normal',
    },

    {
      -logic_name        => 'UpdateMetadata',
      -module            => 'Bio::EnsEMBL::EGPipeline::DNAFeatures::UpdateMetadata',
      -analysis_capacity => 10,
      -parameters        => {},
      -rc_name           => 'normal',
      -flow_into         => WHEN('#email_report#' => ['EmailRepeatReport']),
    },

    {
      -logic_name        => 'EmailRepeatReport',
      -module            => 'Bio::EnsEMBL::EGPipeline::DNAFeatures::EmailRepeatReport',
      -analysis_capacity => 10,
      -parameters        => {
                              email   => $self->o('email'),
                              subject => 'DNA features pipeline: Repeat report for #dbname#',
                              report_dir => $self->o('report_dir'),
                            },
      -rc_name           => 'normal',
    }

  ];
}

1;
