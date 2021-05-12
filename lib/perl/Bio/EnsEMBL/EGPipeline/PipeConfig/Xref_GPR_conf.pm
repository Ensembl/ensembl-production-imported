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

Bio::EnsEMBL::EGPipeline::PipeConfig::Xref_GPR_conf

=head1 DESCRIPTION

Assign Gramene Plant Reactome xrefs available from here:
http://plantreactome.gramene.org/index.php?Itemid=242

=head1 Author

Dan Bolser and James Allen

=cut

package Bio::EnsEMBL::EGPipeline::PipeConfig::Xref_GPR_conf;

use strict;
use warnings;

## EG common configuration (mostly resource classes)
use base ('Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf');

## Hive common configuration (every hive pipeline needs this)
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

use Bio::EnsEMBL::Hive::Version 2.4;

use File::Spec::Functions qw(catdir);

sub default_options {
  my ($self) = @_;
  return {
    %{$self->SUPER::default_options},

    pipeline_name => 'xref_gpr_' . $self->o('ensembl_release'),

    ## We define these here, I think, even though we don't technically
    ## have to to avoid 'first use of x uninitalized' errors during
    ## the pipeline run.
    species      => [],
    division     => [],
    run_all      => 0,
    antispecies  => [],
    meta_filters => {},
    db_type      => 'core',

    ## XRef file provided by Gramene Plant Reactome DB
    xref_reac_file => '',
    xref_path_file => '',

    ## ??
    replace_all => 0,

    external_dbs => {
        'reac' => 'Plant_Reactome_Reaction',
        'path' => 'Plant_Reactome_Pathway',
    },

    ## Not sure if we need the below for this specific pipeline... Are
    ## they are useful because some analysis are optional? Actually I
    ## think we're just being clear about the logic names we'll be
    ## using and their associated modules.
    logic_name => 'gramene_plant_reactome',
    module     => 'Bio::EnsEMBL::EGPipeline::Xref::LoadGPR',

    # use uppercased stable gene ids for mapping
    uppercase_gene_id => 0,

    # Retrieve analysis descriptions from the production database;
    # the supplied registry file will need the relevant server
    # details. Not sure if this is true, because it's defined in
    # Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf
    production_lookup => 1,

    # Entries in the xref table that are not linked to other tables
    # via a foreign key relationship are deleted by default.
    delete_unattached_xref => 1,

    # Default capacity is low, to limit strain on our db servers
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
  ];
}

sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{$self->SUPER::pipeline_wide_parameters},
    'db_type'                => $self->o('db_type'),
    'delete_unattached_xref' => $self->o('delete_unattached_xref'),
    'uppercase_gene_id' => $self->o('uppercase_gene_id'),
  };
}

sub pipeline_analyses {
  my ($self) = @_;

  return [

    {
      -logic_name      => 'SpeciesFactory',
      -module          => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EGSpeciesFactory',
      -input_ids       => [ {} ],
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
                            '2' => ['BackupTables'],
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
                               'object_xref',
                               'xref',
                             ],
                              output_file => catdir($self->o('pipeline_dir'), '#species#', 'pre_pipeline_bkp.sql.gz'),
                            },
      -rc_name           => 'normal',
      -flow_into         => ['SetupGPRAnalysis'],
    },

    {
      -logic_name        => 'SetupGPRAnalysis',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::AnalysisSetup',
      -analysis_capacity => 10,
      -batch_size        => 10,
      -max_retry_count   => 0,
      -parameters        => {
                              logic_name         => $self->o('logic_name'),
                              ## Is the module simply written to the table for book keeping?
                              module             => $self->o('module'),
                              production_lookup  => $self->o('production_lookup'),
                              ## Defined in
                              ## Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf
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
                              'A->1' => ['LoadGPR'],
                            },
    },

    {
      -logic_name      => 'LoadGPR',
      -module          => $self->o('module'),
      -hive_capacity   => $self->o('hive_capacity'),
      -max_retry_count => 0,
      -parameters      => {
                            logic_name        => $self->o('logic_name'),
                            external_dbs      => $self->o('external_dbs'),
                            xref_reac_file    => $self->o('xref_reac_file'),
                            xref_path_file    => $self->o('xref_path_file'),
                            uppercase_gene_id => $self->o('uppercase_gene_id'),
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
 ];
}

1;
