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

Bio::EnsEMBL::EGPipeline::PipeConfig::Xref_VB_conf

=head1 DESCRIPTION

Load VectorBase community-submitted xrefs.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::PipeConfig::Xref_VB_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf');

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use Bio::EnsEMBL::Hive::Version v2.4;
use File::Spec::Functions qw(catdir);

sub default_options {
  my ($self) = @_;
  return {
    %{$self->SUPER::default_options},

    pipeline_name => $self->tagged_pipeline_name('xref_vb', $self->o('ensembl_release')),

    species      => [],
    division     => [],
    run_all      => 0,
    antispecies  => [],
    meta_filters => {},
    db_type      => 'core',
    
    annotations_dir => $self->o('annotation_dir'),

    logic_name => 'brc4_community_annotation',
    module     => 'Bio::EnsEMBL::EGPipeline::Xref::LoadVBCommunityAnnotations',

    vb_external_db       => 'BRC4_Community_Annotation',
    citation_external_db => 'PUBMED',

    # Exclude genes which come from sources with names already
    exclude_name_source => [
      'mirbase_gene',
      'rfam_12.1_gene',
      'trnascan_gene',
    ],
    
    # Exclude genes which come from sources with descriptions already
    exclude_desc_source => [
      'mirbase_gene',
      'rfam_12.1_gene',
      'trnascan_gene',
      'refseq_mdom',
      'refseq_scal',
    ],
    
    description_blacklist => ['Uncharacterized protein', 'AGAP\d.*', 'AAEL\d.*'],
    
    # Retrieve analysis descriptions from the production database;
    # the supplied registry file will need the relevant server details.
    production_lookup => 1,

    # Entries in the xref table that are not linked to other tables
    # via a foreign key relationship are deleted by default.
    delete_unattached_xref => 1,

    # By default, an email is sent for each species when the pipeline
    # is complete, showing the breakdown of xrefs assigned.
    email_xref_report => 1,

    # Default capacity is low, to limit strain on our db servers.
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
    $self->db_cmd("CREATE TABLE gene_descriptions (species varchar(100) NOT NULL, db_name varchar(100) NOT NULL, total int NOT NULL, timing varchar(10))"),
    $self->db_cmd("CREATE TABLE gene_names (species varchar(100) NOT NULL, db_name varchar(100) NOT NULL, total int NOT NULL, timing varchar(10))"),
  ];
}

sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{$self->SUPER::pipeline_wide_parameters},
    'db_type'                => $self->o('db_type'),
    'delete_unattached_xref' => $self->o('delete_unattached_xref'),
    'email_xref_report'      => $self->o('email_xref_report'),
  };
}

sub pipeline_analyses {
  my ($self) = @_;

  return [
    {
      -logic_name        => 'Init',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'mkdir -p ' . $self->o('pipeline_dir'),
      },
      -input_ids  => [{}],
      -rc_name           => 'normal',
      -analysis_capacity => 1,
      -max_retry_count => 0,
      -flow_into  => 'SpeciesFactory',
    },
    {
      -logic_name      => 'SpeciesFactory',
      -module          => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EGSpeciesFactory',
      -parameters      => {
                            species         => $self->o('species'),
                            antispecies     => $self->o('antispecies'),
                            division        => $self->o('division'),
                            run_all         => $self->o('run_all'),
                            meta_filters    => $self->o('meta_filters'),
                            chromosome_flow => 0,
                            regulation_flow => 0,
                            variation_flow  => 0,
                          },
      -max_retry_count => 1,
      -flow_into       => {
                            '2->A' => ['RunPipeline'],
                            'A->2' => ['FinishingTouches'],
                          },
      -rc_name           => 'normal',
    },

    {
      -logic_name      => 'RunPipeline',
      -module          => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -max_retry_count => 0,
      -flow_into       => {
                            '1->A' => WHEN('#email_xref_report#' => ['NamesAndDescriptionsBefore']),
                            'A->1' => ['BackupTables'],
                          },
      -rc_name           => 'normal',
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
      -flow_into         => ['AnalysisSetup']
    },

    {
      -logic_name      => 'AnalysisSetup',
      -module          => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::AnalysisSetup',
      -max_retry_count => 0,
      -parameters      => {
                            logic_name         => $self->o('logic_name'),
                            module             => $self->o('module'),
                            production_lookup  => $self->o('production_lookup'),
                            production_db      => $self->o('production_db'),
                            db_backup_required => 1,
                            db_backup_file     => catdir($self->o('pipeline_dir'), '#species#', 'pre_pipeline_bkp.sql.gz'),
                            delete_existing    => 1,
                            linked_tables      => ['object_xref'],
                          },
      -rc_name           => 'normal',
      -flow_into       => {
                            '1->A' => ['RemoveOrphans'],
                            'A->1' => ['LoadVBCommunityAnnotations'],
                          },
    },

    {
      -logic_name      => 'LoadVBCommunityAnnotations',
      -module          => 'Bio::EnsEMBL::EGPipeline::Xref::LoadVBCommunityAnnotations',
      -max_retry_count => 0,
      -parameters      => {
                            annotation_dir       => $self->o('annotation_dir'),
                            pipeline_dir         => $self->o('pipeline_dir'),
                            logic_name           => $self->o('logic_name'),
                            vb_external_db       => $self->o('vb_external_db'),
                            citation_external_db => $self->o('citation_external_db'),
                            exclude_name_source  => $self->o('exclude_name_source'),
                            exclude_desc_source  => $self->o('exclude_desc_source'),
                          },
      -max_retry_count => 1,
      -hive_capacity   => $self->o('hive_capacity'),
      -rc_name         => 'normal',
    },

    {
      -logic_name      => 'RemoveOrphans',
      -module          => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::SqlCmd',
      -max_retry_count => 0,
      -parameters      => {
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
      -rc_name           => 'normal',
    },

    {
      -logic_name      => 'FinishingTouches',
      -module          => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -max_retry_count => 0,
      -flow_into       => {
                            '1->A' => WHEN('#delete_unattached_xref#' => ['DeleteUnattachedXref']),
                            'A->1' => ['SetupXrefReport'],
                          },
      -rc_name           => 'normal',
    },

    {
      -logic_name      => 'DeleteUnattachedXref',
      -module          => 'Bio::EnsEMBL::EGPipeline::Xref::DeleteUnattachedXref',
      -max_retry_count => 0,
      -parameters      => {},
      -rc_name         => 'normal',
    },

    {
      -logic_name      => 'SetupXrefReport',
      -module          => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -max_retry_count => 0,
      -flow_into       => {
                            '1->A' => WHEN('#email_xref_report#' => ['NamesAndDescriptionsAfter']),
                            'A->1' => WHEN('#email_xref_report#' => ['EmailXrefReport']),
                          },
      -rc_name           => 'normal',
    },

    {
      -logic_name      => 'NamesAndDescriptionsBefore',
      -module          => 'Bio::EnsEMBL::EGPipeline::Xref::NamesAndDescriptions',
      -max_retry_count => 0,
      -parameters      => {
                            timing => 'before',
                          },
      -rc_name         => 'normal',
      -flow_into       => {
                            '2' => ['?table_name=gene_descriptions'],
                            '3' => ['?table_name=gene_names'],
                          }
    },

    {
      -logic_name      => 'NamesAndDescriptionsAfter',
      -module          => 'Bio::EnsEMBL::EGPipeline::Xref::NamesAndDescriptions',
      -max_retry_count => 0,
      -parameters      => {
                            timing => 'after',
                          },
      -rc_name         => 'normal',
      -flow_into       => {
                            '2' => ['?table_name=gene_descriptions'],
                            '3' => ['?table_name=gene_names'],
                          }
    },

    {
      -logic_name      => 'EmailXrefReport',
      -module          => 'Bio::EnsEMBL::EGPipeline::Xref::EmailXrefVBReport',
      -parameters      => {
                            email                => $self->o('email'),
                            subject              => 'Xref (VB) pipeline report for #species#',
                            db_type              => $self->o('db_type'),
                            logic_name           => $self->o('logic_name'),
                            vb_external_db       => $self->o('vb_external_db'),
                            citation_external_db => $self->o('citation_external_db'),
                          },
      -max_retry_count => 1,
      -rc_name         => 'normal',
    },
 ];
}

1;
