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
Bio::EnsEMBL::EGPipeline::PipeConfig::RNAFeatures_conf

=head1 DESCRIPTION

Configuration for running the RNA Features pipeline, which
aligns RNA features against a genome and adds them to a core database.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::PipeConfig::RNAFeatures_conf;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use Bio::EnsEMBL::Hive::Version 2.5;

use base ('Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf');
use File::Spec::Functions qw(catdir);


sub default_options {
  my ($self) = @_;
  
  return {
    %{$self->SUPER::default_options},

    pipeline_name => $self->tagged_pipeline_name('rna_features', $self->o('ensembl_release')),

    # Don't use rel2abs to create the work_dir from $self->o("...") strings.
    # The substitution is done too late. Thus we have duplicated string
    # in the case of the abs-path provided
    
    species => [],
    antispecies => [],
    division => [],
    run_all => 0,
    meta_filters => {},

    # Parameters for dumping and splitting Fasta DNA files...
    max_seq_length          => 10000000,
    max_seq_length_per_file => $self->o('max_seq_length'),
    max_seqs_per_file       => 1000,
    max_files_per_directory => 50,
    max_dirs_per_directory  => $self->o('max_files_per_directory'),

    max_hive_capacity => 50,

    run_cmscan   => 1,
    run_trnascan => 1,
    load_mirbase => 1,

    cmscan_exe => 'cmscan',

    cmscan_cm_file    => {},
    cmscan_logic_name => {},
    cmscan_db_name    => {},
    cmscan_cpu        => 3,
    cmscan_heuristics => 'default',
    cmscan_threshold  => 0.001,
    clanin_file       => undef,
    cmscan_param_hash =>
    {
      cpu            => $self->o('cmscan_cpu'),
      heuristics     => $self->o('cmscan_heuristics'),
      threshold      => $self->o('cmscan_threshold'),
      clanin_file    => $self->o('clanin_file'),
    },
    cmscan_parameters => '',
    cmsscan_resource_class => 'cmscan_4Gb_mem',

    # The blacklist is a foolproof method of excluding Rfam models that
    # you do not want to annotate, perhaps because they generate an excess of
    # alignments. The taxonomic filtering tries to address this, but there may
    # be a small number of models that slip past that filter, but which are
    # nonetheless inappropriate. (An alternative is to ratchet up the strictness
    # of the taxonomic filter, but you may then start excluding appropriate
    # models...)
    rfam_version        => $self->private_conf('RFAM_VERSION'),
    rfam_dir            => catdir($self->private_conf('RFAM_VERSIONS_DIR'), $self->o('rfam_version')),
    rfam_cm_file        => catdir($self->o('rfam_dir'), 'Rfam.cm'),
    rfam_logic_name     => 'cmscan_rfam_'.$self->o('rfam_version'),
    rfam_db_name        => 'RFAM',
    rfam_rrna           => 1,
    rfam_trna           => 0,
    rfam_blacklist      => {
                            'Archaea' =>
                              ['RF00002', 'RF00177', 'RF01118', 'RF01854', 'RF01960', 'RF02514', 'RF02541', 'RF02542', 'RF02543', ],
                            'Bacteria' =>
                              ['RF00002', 'RF00882', 'RF01959', 'RF01960', 'RF02540', 'RF02542', 'RF02543', ],
                            'Eukaryota' =>
                              ['RF00177', 'RF01118', 'RF00169', 'RF01959', 'RF02514', 'RF02540', 'RF02541', 'RF02542', ],
                            'Fungi' =>
                              ['RF00012', 'RF00017', 'RF01848', ],
                            'Metazoa' =>
                              ['RF00882', 'RF00906', 'RF01582', 'RF01675', 'RF01846', 'RF01848', 'RF01849', 'RF01856', 'RF02032', 'RF02625', 'RF02626', 'RF02628', 'RF02647', 'RF02682', ],
                            'Viridiplantae' =>
                              ['RF00012', 'RF00017', 'RF01856', 'RF02628', ],
                            'EnsemblProtists' =>
                              ['RF00012', 'RF00017', 'RF01855', ],
                            'Ensembl' =>
                              ['RF01358', 'RF01376', ],
                            },
    rfam_whitelist      => {
                            'EnsemblProtists' =>
                              ['RF00029', ],
                            },
    rfam_taxonomy_file  => catdir($self->o('rfam_dir'), 'taxonomic_levels.txt'),
    taxonomic_filtering => 1,
    taxonomic_lca       => 0,
    taxonomic_levels    => [],
    taxonomic_threshold => 0.02,
    taxonomic_minimum   => 50,

    # There's not much to choose between cmscan and tRNAscan-SE in terms of
    # annotating tRNA genes, (for some species both produce lots of false
    # positives). If you use tRNAscan-SE, however, you do have the option of
    # including pseudogenes, and you get info about the anticodon in the
    # gene description.
    trnascan_exe        => 'tRNAscan-SE',
    trnascan_logic_name => 'trnascan_align',
    trnascan_db_name    => 'TRNASCAN_SE',
    trnascan_pseudo     => 0,
    trnascan_threshold  => 40,
    trnascan_parameters => '',

    # The annotation from mirBase isn't always available, but if it is,
    # it's useful to load, since it's the definitive source for miRNA.
    mirbase_logic_name => 'mirbase',
    mirbase_db_name    => 'miRBase',
    mirbase_version    => '21',
    mirbase_file       => {},

    analyses =>
    [
      {
        'logic_name'      => $self->o('rfam_logic_name'),
        'db'              => 'Rfam',
        'db_version'      => $self->o('rfam_version'),
        'db_file'         => $self->o('rfam_cm_file'),
        'program'         => 'Infernal',
        'program_version' => '1.1',
        'program_file'    => $self->o('cmscan_exe'),
        'parameters'      => $self->o('cmscan_parameters'),
        'module'          => 'Bio::EnsEMBL::Analysis::Runnable::CMScan',
        'linked_tables'   => ['dna_align_feature'],
      },
      {
        'logic_name'      => $self->o('rfam_logic_name').'_lca',
        'db'              => 'Rfam',
        'db_version'      => $self->o('rfam_version'),
        'db_file'         => $self->o('rfam_cm_file'),
        'program'         => 'Infernal',
        'program_version' => '1.1',
        'program_file'    => $self->o('cmscan_exe'),
        'parameters'      => $self->o('cmscan_parameters'),
        'module'          => 'Bio::EnsEMBL::Analysis::Runnable::CMScan',
        'linked_tables'   => ['dna_align_feature'],
      },
      {
        'logic_name'      => 'cmscan_custom',
        'program'         => 'Infernal',
        'program_version' => '1.1',
        'program_file'    => $self->o('cmscan_exe'),
        'parameters'      => $self->o('cmscan_parameters'),
        'module'          => 'Bio::EnsEMBL::Analysis::Runnable::CMScan',
        'linked_tables'   => ['dna_align_feature'],
      },
      {
        'logic_name'      => $self->o('trnascan_logic_name'),
        'program'         => 'tRNAscan-SE',
        'program_version' => '1.23',
        'program_file'    => $self->o('trnascan_exe'),
        'parameters'      => $self->o('trnascan_parameters'),
        'module'          => 'Bio::EnsEMBL::Analysis::Runnable::tRNAscan',
        'linked_tables'   => ['dna_align_feature'],
      },
      {
        'logic_name'      => $self->o('mirbase_logic_name'),
        'db'              => $self->o('mirbase_db_name'),
        'db_version'      => $self->o('mirbase_version'),
        'module'          => 'Bio::EnsEMBL::EGPipeline::RNAFeatures::miRBase',
        'linked_tables'   => ['dna_align_feature'],
      },
    ],

    # Remove existing DNA features; if => 0 then existing analyses
    # and their features will remain, with the logic_name suffixed by '_bkp'.
    delete_existing => 1,

    # Retrieve analysis descriptions from the production database;
    # the supplied registry file will need the relevant server details.
    production_lookup => 1,

    # An email is sent summarising the alignments for all species,
    # and plots are produced with evalue cut-offs, each with a different colour. 
    evalue_levels =>  {
                        '1e-3' => 'forestgreen',
                        '1e-6' => 'darkorange',
                        '1e-9' => 'firebrick',
                      },
    #
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
    'mkdir -p '.$self->o('pipeline_dir'),
  ];
}

sub pipeline_wide_parameters {
 my ($self) = @_;
 
 return {
   %{$self->SUPER::pipeline_wide_parameters},
   'pipeline_dir' => $self->o('pipeline_dir'),
   'run_cmscan'   => $self->o('run_cmscan'),
   'run_trnascan' => $self->o('run_trnascan'),
   'load_mirbase' => $self->o('load_mirbase'),
   'cmsscan_resource_class' => $self->o('cmsscan_resource_class'),
 };
}

sub pipeline_analyses {
  my ($self) = @_;
  
  my $rfam_logic_name = $self->o('rfam_logic_name');
  if ($self->o('taxonomic_filtering') && $self->o('taxonomic_lca')) {
    $rfam_logic_name .= '_lca';
  }

  return [
    {
      -logic_name        => 'RNAFeatures',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -max_retry_count   => 1,
      -parameters        => {
                            },
      -input_ids         => [ {} ],
      -flow_into         => {
                              '1->A' => ['dbFactory_features'],
                              'A->1' => ['DbFactory_Report'],
                            },
      -meadow_type       => 'LOCAL',
    },

    {
      -logic_name      => 'dbFactory_features',
      -module          => 'Bio::EnsEMBL::Production::Pipeline::Common::DbFactory',
      -max_retry_count => 0,
      -parameters      => {
                            species      => $self->o('species'),
                            division     => $self->o('division'),
                            antispecies  => $self->o('antispecies'),
                            run_all      => $self->o('run_all'),
                            meta_filters => $self->o('meta_filters'),
                          },
      -rc_name         => 'default',
      -flow_into       => {
                            '2' => ['CreateWD'],
                          },
    },

    {
      -logic_name        => 'CreateWD',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters         => {
                               cmscan_wd  => catdir($self->o('pipeline_dir'), '#species#', 'cmscan_wd'),
                               trnascan_wd => catdir($self->o('pipeline_dir'), '#species#', 'trnascan_wd'),
                               cmd => 'mkdir -p #cmscan_wd# #trnascan_wd#',
      },
      -rc_name           => 'normal',
      -max_retry_count => 0,
      -flow_into         => ['BackupDatabase'],
    },

    {
      -logic_name        => 'BackupDatabase',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::DatabaseDumper',
      -analysis_capacity => 5,
      -max_retry_count   => 1,
      -parameters        => {
                              table_list  => [
                                'analysis',
                                'analysis_description',
                                'dna_align_feature',
                                'dna_align_feature_attrib',
                              ],
                              output_file => catdir($self->o('pipeline_dir'), '#species#', 'pre_pipeline_bkp.sql.gz'),
                            },
      -rc_name           => 'normal',
      -flow_into         => {
                              '1->A' => ['RNAAnalysisFactory'],
                              'A->1' => ['DbAwareSpeciesFactoryFeatures'],
                            },
    },
    {
      -logic_name        => 'RNAAnalysisFactory',
      -module            => 'Bio::EnsEMBL::EGPipeline::RNAFeatures::RNAAnalysisFactory',
      -max_retry_count   => 0,
      -batch_size        => 10,
      -parameters        => {
                              analyses            => $self->o('analyses'),
                              rfam_logic_name     => $rfam_logic_name,
                              cmscan_logic_name   => $self->o('cmscan_logic_name'),
                              cmscan_cm_file      => $self->o('cmscan_cm_file'),
                              trnascan_logic_name => $self->o('trnascan_logic_name'),
                              mirbase_logic_name  => $self->o('mirbase_logic_name'),
                              mirbase_files       => $self->o('mirbase_file'),
                              pipeline_dir        => $self->o('pipeline_dir'),
                              db_backup_file      => catdir($self->o('pipeline_dir'), '#species#', 'pre_pipeline_bkp.sql.gz'),
                            },
      -rc_name           => 'normal',
      -flow_into         => {
                              '2' => ['AnalysisSetup'],
                            },
    },
    {
      -logic_name        => 'DbAwareSpeciesFactoryFeatures',
      -module            => 'Bio::EnsEMBL::Production::Pipeline::Common::DbAwareSpeciesFactory',
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
      -rc_name           => 'normal',
      -flow_into         => {
                              '2->A' => WHEN('#run_cmscan#' => ['TaxonomicFilter']),
                              'A->2' => ['dbAwareMetaCoordDummy'],
                            },
    },
    {
      -logic_name        => 'dbAwareMetaCoordDummy',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -max_retry_count   => 1,
      -can_be_empty      => 1,
      -parameters        => {},
      -rc_name           => 'normal',
      -flow_into         => {
                              '1->A' => ['AnnotateRNAFeatures'],
                              'A->1' => ['MetaCoords'],
                            },
      -meadow_type       => 'LOCAL',
    },

    {
      -logic_name        => 'MetaCoords',
      -module            => 'Bio::EnsEMBL::EGPipeline::CoreStatistics::MetaCoords',
      -analysis_capacity => 10,
      -batch_size        => 10,
      -max_retry_count   => 1,
      -parameters        => {},
      -rc_name           => 'normal',
    },
    {
      -logic_name        => 'AnalysisSetup',
      -module            => 'Bio::EnsEMBL::Production::Pipeline::Common::AnalysisSetup', 
      -analysis_capacity => 10,
      -batch_size        => 10,
      -max_retry_count   => 0,
      -parameters        => {
                              db_backup_required => 1,
                              delete_existing    => $self->o('delete_existing'),
                              production_lookup  => $self->o('production_lookup'),
                              production_db      => $self->o('production_db'),
                            },
      -meadow_type       => 'LOCAL',
      -flow_into         => ['DeleteAttributes'],
    },

    {
      -logic_name      => 'DeleteAttributes',
      -module          => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::SqlCmd',
      -max_retry_count => 0,
      -parameters      => {
                             sql => [
                               'DELETE dafa.* FROM '.
                                 'dna_align_feature_attrib dafa LEFT OUTER JOIN '.
                                 'dna_align_feature daf USING (dna_align_feature_id) '.
                                 'WHERE daf.dna_align_feature_id IS NULL',
                             ]
                           },
      -meadow_type     => 'LOCAL',
    },

    {
      -logic_name        => 'TaxonomicFilter',
      -module            => 'Bio::EnsEMBL::EGPipeline::RNAFeatures::TaxonomicFilter',
      -analysis_capacity => 100,
      -batch_size        => 100,
      -max_retry_count   => 3,
      -parameters        => {
                              rfam_cm_file        => $self->o('rfam_cm_file'),
                              filtered_cm_file    => catdir($self->o('pipeline_dir'), '#species#', 'Rfam.filtered.cm'),
                              rfam_logic_name     => $rfam_logic_name,
                              rfam_rrna           => $self->o('rfam_rrna'),
                              rfam_trna           => $self->o('rfam_trna'),
                              rfam_blacklist      => $self->o('rfam_blacklist'),
                              rfam_whitelist      => $self->o('rfam_whitelist'),
                              rfam_taxonomy_file  => $self->o('rfam_taxonomy_file'),
                              taxonomic_filtering => $self->o('taxonomic_filtering'),
                              taxonomic_lca       => $self->o('taxonomic_lca'),
                              taxonomic_levels    => $self->o('taxonomic_levels'),
                              taxonomic_threshold => $self->o('taxonomic_threshold'),
                              taxonomic_minimum   => $self->o('taxonomic_minimum'),
                            },
      -rc_name           => '4Gb_mem',
      -flow_into         => ['CMScanIndex'],
    },

    {
      -logic_name        => 'CMScanIndex',
      -module            => 'Bio::EnsEMBL::EGPipeline::RNAFeatures::CMScanIndex',
      -analysis_capacity => 30,
      -batch_size        => 100,
      -max_retry_count   => 1,
      -parameters        => {
                              rfam_cm_file      => catdir($self->o('pipeline_dir'), '#species#', 'Rfam.filtered.cm'),
                              rfam_logic_name   => $rfam_logic_name,
                              cmscan_cm_file    => $self->o('cmscan_cm_file'),
                              cmscan_logic_name => $self->o('cmscan_logic_name'),
                              parameters_hash   => $self->o('cmscan_param_hash'),
                            },
      -rc_name           => 'normal',
    },

    {
      -logic_name        => 'AnnotateRNAFeatures',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -max_retry_count   => 1,
      -parameters        => {},
      -flow_into         => WHEN(
                              '#run_cmscan# or #run_trnascan#' => ['DumpGenome'],
                              '#load_mirbase#'                 => ['miRBase'],
                            ),
      -meadow_type       => 'LOCAL',
    },

    {
      -logic_name        => 'DumpGenome',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::DumpGenome',
      -analysis_capacity => 100,
      -batch_size        => 100,
      -parameters        => {
                              genome_dir => catdir($self->o('pipeline_dir'), '#species#'),
                            },
      -rc_name           => 'normal',
      -flow_into         => ['SplitDumpFile'],
    },

    {
      -logic_name        => 'SplitDumpFile',
      -module            => 'Bio::EnsEMBL::EGPipeline::RNAFeatures::SplitDumpFile',
      -analysis_capacity => 50,
      -batch_size        => 50,
      -max_retry_count   => 1,
      -parameters        => {
                              fasta_file              => '#genome_file#',
                              max_seq_length_per_file => $self->o('max_seq_length_per_file'),
                              max_seqs_per_file       => $self->o('max_seqs_per_file'),
                              max_files_per_directory => $self->o('max_files_per_directory'),
                              max_dirs_per_directory  => $self->o('max_dirs_per_directory'),
                            },
      -rc_name           => '8Gb_mem',
      -flow_into         => {
                              '3' => ['CMScanFactory'],
                              '4' => ['tRNAscan'],
                            },
    },

    {
      -logic_name        => 'CMScanFactory',
      -module            => 'Bio::EnsEMBL::EGPipeline::RNAFeatures::CMScanFactory',
      -analysis_capacity => 50,
      -batch_size        => 100,
      -max_retry_count   => 1,
      -parameters        => {
                              rfam_cm_file      => catdir($self->o('pipeline_dir'), '#species#', 'Rfam.filtered.cm'),
                              rfam_logic_name   => $rfam_logic_name,
                              cmscan_cm_file    => $self->o('cmscan_cm_file'),
                              cmscan_logic_name => $self->o('cmscan_logic_name'),
                              cmscan_db_name    => $self->o('cmscan_db_name'),
                              parameters_hash   => $self->o('cmscan_param_hash'),
                              max_seq_length    => $self->o('max_seq_length'),
                            },
      -rc_name           => '8Gb_mem',
      -flow_into         => {
                              '2' => ['CMScan'],
                            },
    },

    {
      -logic_name        => 'CMScan',
      -module            => 'Bio::EnsEMBL::EGPipeline::RNAFeatures::CMScan',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -max_retry_count   => 1,
      -analysis_capacity => 50,
      -batch_size        => 10,
      -parameters        => {
                              escape_branch => -1,
                              db_name       => $self->o('rfam_db_name'),
                              workdir       => catdir($self->o('pipeline_dir'), '#species#', 'cmscan_wd'),
                            },
      -rc_name           => $self->o('cmsscan_resource_class'),
      -flow_into         => {
                              '-1' => ['CMScan_HighMem'],
                            },
    },

    {
      -logic_name        => 'CMScan_HighMem',
      -module            => 'Bio::EnsEMBL::EGPipeline::RNAFeatures::CMScan',
      -can_be_empty      => 1,
      -hive_capacity     => $self->o('max_hive_capacity'),
      -max_retry_count   => 0,
      -parameters        => {
                              db_name => $self->o('rfam_db_name'),
                              workdir       => catdir($self->o('pipeline_dir'), '#species#', 'cmscan_wd'),
                            },
      -rc_name           => 'cmscan_8Gb_mem',
      -flow_into         => {
                              '-1' => ['CMScan_HighMem_16'],
                            },
    },

    {
      -logic_name        => 'CMScan_HighMem_16',
      -module            => 'Bio::EnsEMBL::EGPipeline::RNAFeatures::CMScan',
      -can_be_empty      => 1,
      -hive_capacity     => $self->o('max_hive_capacity'),
      -max_retry_count   => 1,
      -parameters        => {
                              db_name => $self->o('rfam_db_name'),
                              workdir       => catdir($self->o('pipeline_dir'), '#species#', 'cmscan_wd'),
                            },
      -rc_name           => 'cmscan_16Gb_mem',
    },

    {
      -logic_name        => 'tRNAscan',
      -module            => 'Bio::EnsEMBL::EGPipeline::RNAFeatures::tRNAscan',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -analysis_capacity => 20,
      -batch_size        => 50,
      -max_retry_count   => 3,
      -parameters        => {
                              logic_name   => $self->o('trnascan_logic_name'),
                              db_name      => $self->o('trnascan_db_name'),
                              pseudo       => $self->o('trnascan_pseudo'),
                              threshold    => $self->o('trnascan_threshold'),
                              workdir      => catdir($self->o('pipeline_dir'), '#species#', 'trnascan_wd'),
                            },
      -rc_name           => 'normal',
    },

    {
      -logic_name        => 'miRBase',
      -module            => 'Bio::EnsEMBL::EGPipeline::RNAFeatures::miRBase',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -max_retry_count   => 1,
      -parameters        => {
                              logic_name => $self->o('mirbase_logic_name'),
                              db_name    => $self->o('mirbase_db_name'),
                              files      => $self->o('mirbase_file'),
                            },
      -rc_name           => 'normal',
    },
    {
      -logic_name        => 'DbFactory_Report',
      -module            => 'Bio::EnsEMBL::Production::Pipeline::Common::DbFactory',
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
      -flow_into         => {
                              '2->A' => ['FetchAlignments'],
                              'A->1' => ['SummariseAlignments'],
                            },
      -meadow_type       => 'LOCAL',
    },
    {
      -logic_name        => 'FetchAlignments',
      -module            => 'Bio::EnsEMBL::EGPipeline::RNAFeatures::FetchAlignments',
      -analysis_capacity => 10,
      -batch_size        => 100,
      -max_retry_count   => 1,
      -parameters        => {
                              rfam_logic_name     => $rfam_logic_name,
                              cmscan_cm_file      => $self->o('cmscan_cm_file'),
                              cmscan_logic_name   => $self->o('cmscan_logic_name'),
                              trnascan_logic_name => $self->o('trnascan_logic_name'),
                              mirbase_logic_name  => $self->o('mirbase_logic_name'),
                              alignment_dir       => catdir($self->o('pipeline_dir'), '#species#'),
                            },
      -rc_name           => 'normal',
    },

    {
      -logic_name        => 'SummariseAlignments',
      -module            => 'Bio::EnsEMBL::EGPipeline::RNAFeatures::SummariseAlignments',
      -analysis_capacity => 10,
      -batch_size        => 10,
      -max_retry_count   => 1,
      -parameters        => {
                              pipeline_dir  => $self->o('pipeline_dir'),
                              evalue_levels => $self->o('evalue_levels'),
                            },
      -flow_into         => ['EmailRNAReport'],
      -meadow_type       => 'LOCAL',
    },

    {
      -logic_name        => 'EmailRNAReport',
      -module            => 'Bio::EnsEMBL::EGPipeline::RNAFeatures::EmailRNAReport',
      -analysis_capacity => 10,
      -batch_size        => 10,
      -max_retry_count   => 1,
      -parameters        => {
                              email              => $self->o('email'),
                              subject            => 'RNA features pipeline report',
                              cmscan_threshold   => $self->o('cmscan_threshold'),
                              trnascan_threshold => $self->o('trnascan_threshold'),
                              pipeline_dir       => $self->o('pipeline_dir'),
                              evalue_levels      => $self->o('evalue_levels'),
                            },
      -rc_name           => 'normal',
    },

  ];
}

sub resource_classes {
  my ($self) = @_;

  my $queue = $self->o('queue_name');
  my @mems = (4, 8, 16, 32);
  my $cpu = $self->o('cmscan_cpu');
  my $time = "24:00:00";

  my %resources = %{$self->SUPER::resource_classes};

  for my $mem (@mems) {
    my $name = "cmscan_${mem}Gb_mem";
    $resources{$name} = $self->make_resource({
      queue => $queue,
      memory => $mem * 1000,
      temp_memory => $mem * 1000,
      cpus => $cpu,
      time => $time
    });
  }
  return \%resources;
}

1;
