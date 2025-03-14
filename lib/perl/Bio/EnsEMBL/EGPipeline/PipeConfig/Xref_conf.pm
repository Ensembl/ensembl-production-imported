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

Bio::EnsEMBL::EGPipeline::PipeConfig::Xref_conf

=head1 DESCRIPTION

Assign UniParc and UniParc-derived xrefs.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::PipeConfig::Xref_conf;

use strict;
use warnings;

## EG common configuration (mostly resource classes)
use base ('Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf');

## Hive common configuration (every hive pipeline needs this)
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

use Bio::EnsEMBL::Hive::Version v2.4;

use File::Spec::Functions qw(catdir);

sub default_options {
  my ($self) = @_;
  return {
    %{$self->SUPER::default_options},

    pipeline_name => 'xref_' . $self->o('ensembl_release'),

    species      => [],
    division     => [],
    run_all      => 0,
    antispecies  => [],
    meta_filters => {},
    db_type      => 'core',

    # UniParc and Uniprot production databases are for internal use only
    #  and can not be accessed externally. Should be replaced with REST API calls.
    remote_uniparc_db => $self->private_conf('ENSEMBL_REMOTE_UNIPARC_DB'),
    remote_uniprot_db => $self->private_conf('ENSEMBL_REMOTE_UNIPPROT_DB'),
    uniparc_dbm_cache_dir => $self->private_conf('UNIPARC_DBM_CACHE_DIR'),
    uniparc_dbm_cache_name =>  'uniparc_cache.tkh',

    replace_all           => 0,
    gene_name_source      => [],
    overwrite_gene_name   => 0,
    description_source    => [],
    overwrite_description => 0,
    description_blacklist => ['Uncharacterized protein', 'AGAP\d.*', 'AAEL\d.*'],

    load_uniprot        => 1,
    load_uniprot_go     => 1,
    load_uniprot_xrefs  => 1,

    uniparc_external_db   => 'UniParc',
    uniprot_external_dbs  => {
      'reviewed'   => 'Uniprot/SWISSPROT',
      'unreviewed' => 'Uniprot/SPTREMBL',
      'splicevar'  => 'Uniprot/Varsplic',
    },
    uniprot_gn_external_db => 'Uniprot_gn',
    uniprot_go_external_db => 'GO',
    uniprot_xref_external_dbs => {
      'ChEMBL'          => 'ChEMBL',
      'EMBL'            => 'EMBL',
      'GeneID'          => 'EntrezGene',
      'MEROPS'          => 'MEROPS',
      'PDB'             => 'PDB',
      'RefSeq'          => 'RefSeq_peptide',
      'STRING'          => 'STRING',
    },
    uniprot_secondary_external_dbs => {
      'RefSeq_peptide' => 'RefSeq_dna',
    },

    checksum_logic_name => 'xrefchecksum',
    checksum_module     => 'Bio::EnsEMBL::EGPipeline::Xref::LoadUniParc',

    uniparc_transitive_logic_name => 'xrefuniparc',
    uniparc_transitive_module     => 'Bio::EnsEMBL::EGPipeline::Xref::LoadUniProt',

    uniprot_transitive_logic_name  => 'xrefuniprot',
    uniprot_transitive_xref_module => 'Bio::EnsEMBL::EGPipeline::Xref::LoadUniProtXrefs',
    
    uniprot_transitive_go_logic_name  => 'gouniprot',
    uniprot_transitive_go_module      => 'Bio::EnsEMBL::EGPipeline::Xref::LoadUniProtGO',
    
    # Retrieve analysis descriptions from the production database;
    # the supplied registry file will need the relevant server details.
    production_lookup => 1,

    # Entries in the xref table that are not linked to other tables
    # via a foreign key relationship are deleted by default.
    delete_unattached_xref => 1,

    # By default, an email is sent for each species when the pipeline
    # is complete, showing the breakdown of xrefs assigned.
    email_xref_report => 1,
    
    # Default capacity is low, to limit strain on our db servers and
    # UniProt's.
    hive_capacity => 10,
  }
}

sub beekeeper_extra_cmdline_options {
  my ($self) = @_;

  my $options = join(' ',
    $self->SUPER::beekeeper_extra_cmdline_options,
    "-reg_conf ".$self->o('registry')
  );

  return $options;
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
    'mkdir -p '.$self->o('pipeline_dir'),
    $self->db_cmd("CREATE TABLE gene_descriptions (species varchar(100) NOT NULL, db_name varchar(100) NOT NULL, total int NOT NULL, timing varchar(10))"),
    $self->db_cmd("CREATE TABLE gene_names (species varchar(100) NOT NULL, db_name varchar(100) NOT NULL, total int NOT NULL, timing varchar(10))"),
  ];
}

sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{$self->SUPER::pipeline_wide_parameters},
    'db_type'                => $self->o('db_type'),
    'load_uniprot'           => $self->o('load_uniprot'),
    'load_uniprot_go'        => $self->o('load_uniprot_go'),
    'load_uniprot_xrefs'     => $self->o('load_uniprot_xrefs'),
    'delete_unattached_xref' => $self->o('delete_unattached_xref'),
    'email_xref_report'      => $self->o('email_xref_report'),
  };
}

sub pipeline_analyses {
  my ($self) = @_;

  return [
    {
      -logic_name      => 'InitialiseXref',
      -module          => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -input_ids       => [ {} ],
      -max_retry_count => 0,
      -flow_into       => {
                            '1->A' => ['ImportUniParc'],
                            'A->1' => ['SpeciesFactory'],
                          },
      -meadow_type     => 'LOCAL',
    },
    
    {
      -logic_name      => 'ImportUniParc',
      -module          => 'Bio::EnsEMBL::EGPipeline::Xref::ImportUniParc',
      -max_retry_count => 1,
      -parameters      => {
                            uniparc_dbm_cache_dir => $self->o('uniparc_dbm_cache_dir'),
                            uniparc_dbm_cache_name => $self->o('uniparc_dbm_cache_name'),
                            dbm_create_script => catdir($self->o('ensembl_production_imported_scripts_dir'), 'uniparc_index', 'create_uniparc_dbm.py'),
                          },
      -rc_name         => 'datamove_64Gb_mem',
    },

    {
      -logic_name      => 'SpeciesFactory',
      -module          => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EGSpeciesFactory',
      -max_retry_count => 1,
      -parameters      => {
                            species         => $self->o('species'),
                            antispecies     => $self->o('antispecies'),
                            division        => $self->o('division'),
                            run_all         => $self->o('run_all'),
                            meta_filters    => $self->o('meta_filters'),

                            ## These are used by the species factory
                            ## to flow into specific analyses for
                            ## specific species data. i.e. if a
                            ## species has chromosomes or not, you may
                            ## want to performe a specific analysis.
                            chromosome_flow => 0,
                            regulation_flow => 0,
                            variation_flow  => 0,
                          },
      -flow_into       => {
                            '2->A' => ['RunPipeline'],
                            'A->2' => WHEN('#email_xref_report#' => ['SetupXrefReport']),
                          },
      -meadow_type     => 'LOCAL',
    },

    {
      -logic_name      => 'RunPipeline',
      -module          => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -max_retry_count => 0,
      -flow_into       => {
                            '1->A' => WHEN('#email_xref_report#' => ['NamesAndDescriptionsBefore']),
                            'A->1' => ['BackupTables'],
                          },
      -meadow_type     => 'LOCAL',
    },

    {
      -logic_name        => 'BackupTables',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::DatabaseDumper',
      -analysis_capacity => 5,
      -max_retry_count   => 1,
      -parameters        => {
                             table_list => [
                               'analysis',
                               'analysis_description',
                               'dependent_xref',
                               'gene',
                               'identity_xref',
                               'interpro',
                               'object_xref',
                               'ontology_xref',
                               'xref',
                             ],
                              output_file => catdir($self->o('pipeline_dir'), '#species#', 'pre_pipeline_bkp.sql.gz'),
                            },
      -rc_name           => 'normal',
      -flow_into         => ['SetupUniParc']
    },

    {
      -logic_name        => 'SetupUniParc',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::AnalysisSetup',
      -analysis_capacity => 10,
      -batch_size        => 10,
      -max_retry_count   => 0,
      -parameters        => {
                              logic_name         => $self->o('checksum_logic_name'),
                              module             => $self->o('checksum_module'),
                              production_lookup  => $self->o('production_lookup'),
                              production_db      => $self->o('production_db'),
                              db_backup_required => 1,
                              db_backup_file     => catdir($self->o('pipeline_dir'), '#species#', 'pre_pipeline_bkp.sql.gz'),
                              delete_existing    => 1,
                              linked_tables      => ['object_xref'],
                              output_logic_name  => 1,
                            },
      -meadow_type       => 'LOCAL',
      -flow_into         => {
                              '1->A' => ['RemoveOrphans'],
                              'A->1' => ['LoadUniParc'],
                            },
    },

    {
      -logic_name      => 'LoadUniParc',
      -module          => $self->o('checksum_module'),
      -hive_capacity   => $self->o('hive_capacity'),
      -max_retry_count => 0,
      -parameters      => {
                            uniparc_dbm_cache => catdir($self->o('uniparc_dbm_cache_dir'), $self->o('uniparc_dbm_cache_name')),
                            dbm_query_script => catdir($self->o('ensembl_production_imported_scripts_dir'), 'uniparc_index', 'query_uniparc_dbm.py'),
                            upi_query_dir      => catdir($self->o('pipeline_dir'), '#species#', 'upi_query'),
                            logic_name  => $self->o('checksum_logic_name'),
                            external_db => $self->o('uniparc_external_db'),
                          },
      -rc_name         => 'normal',
      -flow_into       => WHEN('#load_uniprot#' => ['SetupUniProt']),
    },

    {
      -logic_name        => 'SetupUniProt',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::AnalysisSetup',
      -analysis_capacity => 10,
      -batch_size        => 10,
      -max_retry_count   => 0,
      -parameters        => {
                              logic_name         => $self->o('uniparc_transitive_logic_name'),
                              module             => $self->o('uniparc_transitive_module'),
                              production_lookup  => $self->o('production_lookup'),
                              production_db      => $self->o('production_db'),
                              db_backup_required => 1,
                              db_backup_file     => catdir($self->o('pipeline_dir'), '#species#', 'pre_pipeline_bkp.sql.gz'),
                              delete_existing    => 1,
                              linked_tables      => ['object_xref'],
                              output_logic_name  => 1,
                            },
      -meadow_type       => 'LOCAL',
      -flow_into         => {
                              '1->A' => ['RemoveOrphans'],
                              'A->1' => ['LoadUniProt'],
                            },
    },

    {
      -logic_name      => 'LoadUniProt',
      -module          => $self->o('uniparc_transitive_module'),
      -hive_capacity   => $self->o('hive_capacity'),
      -max_retry_count => 0,
      -parameters      => {
                            uniparc_db            => $self->o('remote_uniparc_db'),
                            uniprot_db            => $self->o('remote_uniprot_db'),
                            replace_all           => $self->o('replace_all'),
                            gene_name_source      => $self->o('gene_name_source'),
                            overwrite_gene_name   => $self->o('overwrite_gene_name'),
                            description_source    => $self->o('description_source'),
                            overwrite_description => $self->o('overwrite_description'),
                            description_blacklist => $self->o('description_blacklist'),
                            logic_name            => $self->o('uniparc_transitive_logic_name'),
                            external_dbs          => $self->o('uniprot_external_dbs'),
                          },
      -rc_name         => 'normal',
      -flow_into       => [
                            WHEN('#load_uniprot_go#'    => ['SetupUniProtGO']),
                            WHEN('#load_uniprot_xrefs#' => ['SetupUniProtXrefs']),
                          ],
    },

    {
      -logic_name        => 'SetupUniProtGO',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::AnalysisSetup',
      -analysis_capacity => 10,
      -batch_size        => 10,
      -max_retry_count   => 0,
      -parameters        => {
                              logic_name         => $self->o('uniprot_transitive_go_logic_name'),
                              module             => $self->o('uniprot_transitive_go_module'),
                              production_lookup  => $self->o('production_lookup'),
                              production_db      => $self->o('production_db'),
                              db_backup_required => 1,
                              db_backup_file     => catdir($self->o('pipeline_dir'), '#species#', 'pre_pipeline_bkp.sql.gz'),
                              delete_existing    => 1,
                              linked_tables      => ['object_xref'],
                              output_logic_name  => 1,
                            },
      -meadow_type       => 'LOCAL',
      -flow_into         => {
                              '1->A' => ['RemoveOrphans'],
                              'A->1' => ['LoadUniProtGO'],
                            },
    },

    {
      -logic_name      => 'LoadUniProtGO',
      -module          => $self->o('uniprot_transitive_go_module'),
      -hive_capacity   => $self->o('hive_capacity'),
      -max_retry_count => 0,
      -parameters      => {
                            uniprot_db           => $self->o('remote_uniprot_db'),
                            replace_all          => $self->o('replace_all'),
                            logic_name           => $self->o('uniprot_transitive_go_logic_name'),
                            external_db          => $self->o('uniprot_go_external_db'),
                            uniprot_external_dbs => $self->o('uniprot_external_dbs'),
                          },
      -rc_name         => 'normal',
    },

    {
      -logic_name        => 'SetupUniProtXrefs',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::AnalysisSetup',
      -analysis_capacity => 10,
      -batch_size        => 10,
      -max_retry_count   => 0,
      -parameters        => {
                              logic_name         => $self->o('uniprot_transitive_logic_name'),
                              module             => $self->o('uniprot_transitive_xref_module'),
                              production_lookup  => $self->o('production_lookup'),
                              production_db      => $self->o('production_db'),
                              db_backup_required => 1,
                              db_backup_file     => catdir($self->o('pipeline_dir'), '#species#', 'pre_pipeline_bkp.sql.gz'),
                              delete_existing    => 1,
                              linked_tables      => ['object_xref'],
                              output_logic_name  => 1,
                            },
      -meadow_type       => 'LOCAL',
      -flow_into         => {
                              '1->A' => ['RemoveOrphans'],
                              'A->1' => ['LoadUniProtXrefs'],
                            },
    },

    {
      -logic_name      => 'LoadUniProtXrefs',
      -module          => $self->o('uniprot_transitive_xref_module'),
      -hive_capacity   => $self->o('hive_capacity'),
      -max_retry_count => 0,
      -parameters      => {
                            uniprot_db             => $self->o('remote_uniprot_db'),
                            replace_all            => $self->o('replace_all'),
                            logic_name             => $self->o('uniprot_transitive_logic_name'),
                            external_dbs           => $self->o('uniprot_xref_external_dbs'),
                            secondary_external_dbs => $self->o('uniprot_secondary_external_dbs'),
                            uniprot_external_dbs   => $self->o('uniprot_external_dbs'),
                          },
      -rc_name         => 'normal',
    },

    {
      -logic_name        => 'RemoveOrphans',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::SqlCmd',
      -analysis_capacity => 10,
      -batch_size        => 10,
      -max_retry_count   => 0,
      -parameters        => {
                               sql => [
                                 'DELETE dx.* FROM '.
                                   'dependent_xref dx LEFT OUTER JOIN '.
                                   'object_xref ox USING (object_xref_id) '.
                                   'WHERE ox.object_xref_id IS NULL',
                                 'DELETE onx.* FROM '.
                                   'ontology_xref onx LEFT OUTER JOIN '.
                                   'object_xref ox USING (object_xref_id) '.
                                   'WHERE ox.object_xref_id IS NULL',
                                 'UPDATE gene g LEFT OUTER JOIN '.
                                   'xref x ON g.display_xref_id = x.xref_id '.
                                   'SET g.display_xref_id = NULL '.
                                   'WHERE x.xref_id IS NULL',
                               ]
                             },
      -meadow_type       => 'LOCAL',
      -flow_into         => WHEN('#delete_unattached_xref#' => ['DeleteUnattachedXref']),
    },

    {
      -logic_name        => 'DeleteUnattachedXref',
      -module            => 'Bio::EnsEMBL::EGPipeline::Xref::DeleteUnattachedXref',
      -analysis_capacity => 10,
      -batch_size        => 10,
      -max_retry_count   => 0,
      -parameters        => {},
      -rc_name           => 'normal',
    },

    {
      -logic_name        => 'SetupXrefReport',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -analysis_capacity => 10,
      -batch_size        => 10,
      -max_retry_count   => 0,
      -flow_into         => {
                              '1->A' => ['NamesAndDescriptionsAfter'],
                              'A->1' => ['EmailXrefReport'],
                            },
      -meadow_type       => 'LOCAL',
    },

    {
      -logic_name        => 'NamesAndDescriptionsBefore',
      -module            => 'Bio::EnsEMBL::EGPipeline::Xref::NamesAndDescriptions',
      -analysis_capacity => 10,
      -batch_size        => 10,
      -max_retry_count   => 0,
      -parameters        => {
                              timing => 'before',
                            },
      -rc_name           => 'normal',
      -flow_into         => {
                              '2' => ['?table_name=gene_descriptions'],
                              '3' => ['?table_name=gene_names'],
                            }
    },

    {
      -logic_name        => 'NamesAndDescriptionsAfter',
      -module            => 'Bio::EnsEMBL::EGPipeline::Xref::NamesAndDescriptions',
      -analysis_capacity => 10,
      -batch_size        => 10,
      -max_retry_count   => 0,
      -parameters        => {
                              timing => 'after',
                            },
      -rc_name           => 'normal',
      -flow_into         => {
                              '2' => ['?table_name=gene_descriptions'],
                              '3' => ['?table_name=gene_names'],
                            }
    },

    {
      -logic_name        => 'EmailXrefReport',
      -module            => 'Bio::EnsEMBL::EGPipeline::Xref::EmailXrefReport',
      -analysis_capacity => 10,
      -batch_size        => 10,
      -max_retry_count   => 1,
      -parameters        => {
                              email                            => $self->o('email'),
                              subject                          => 'Xref pipeline report for #species#',
                              db_type                          => $self->o('db_type'),
                              load_uniprot                     => $self->o('load_uniprot'),
                              load_uniprot_go                  => $self->o('load_uniprot_go'),
                              load_uniprot_xrefs               => $self->o('load_uniprot_xrefs'),
                              checksum_logic_name              => $self->o('checksum_logic_name'),
                              uniparc_transitive_logic_name    => $self->o('uniparc_transitive_logic_name'),
                              uniprot_transitive_logic_name    => $self->o('uniprot_transitive_logic_name'),
                              uniprot_transitive_go_logic_name => $self->o('uniprot_transitive_go_logic_name'),
                              uniparc_external_db              => $self->o('uniparc_external_db'),
                              uniprot_external_dbs             => $self->o('uniprot_external_dbs'),
                              uniprot_go_external_db           => $self->o('uniprot_go_external_db'),
                              uniprot_xref_external_dbs        => $self->o('uniprot_xref_external_dbs'),
                              replace_all                      => $self->o('replace_all'),
                              gene_name_source                 => $self->o('gene_name_source'),
                              overwrite_gene_name              => $self->o('overwrite_gene_name'),
                              description_source               => $self->o('description_source'),
                              overwrite_description            => $self->o('overwrite_description'),
                            },
      -rc_name           => 'normal',
    },
 ];
}

1;
