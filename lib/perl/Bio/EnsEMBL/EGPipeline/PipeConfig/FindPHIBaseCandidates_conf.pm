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

  Bio::EnsEMBL::EGPipeline::PipeConfig::FindPHIBaseCandidates_conf

=head1 SYNOPSIS


  init_pipeline.pl Bio::EnsEMBL::EGPipeline::PipeConfig::FindPHIBaseCandidates_conf  -pipeline_url $EHIVE_URL  -registry $REGISTRY_FILE -blast_db_dir $BLAST_DB_DIRECTORY -input_file $INPUT_ROTHAMSTED_FILE -hive_force_init 1

  runWorker.pl -url $EHIVE_URL

  runWorker.pl -url $EHIVE_URL --reg_conf $REGISTRY_FILE



=head1 DESCRIPTION

  This is an example pipeline put together from five basic building blocks:

  Analysis_1: JobFactory.pm is used to turn the list of files in a given directory into jobs

      these jobs are sent down the branch #2 into the second analysis

  Analysis_2: JobFactory.pm is used to run a wc command to determine the size of a file
              (and format the output with sed), then capture the command's object, putting
              the file size into a parameter for later use.

  Analysis_3: SystemCmd.pm is used to run these compression/decompression jobs in parallel.

  Analysis_4: JobFactory.pm is used to run a wc command to determine the size of a file
              (and format the output with sed), then capture the command's object, putting
              the file size into a parameter for later use.

  Analysis_5: SystemCmd.pm is used to run the notify-send command, displaying a message on the screen.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <https://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <https://www.ensembl.org/Help/Contact>.

=cut

package Bio::EnsEMBL::EGPipeline::PipeConfig::FindPHIBaseCandidates_conf;


use strict;
use warnings;

## EG common configuration (mostly resource classes)
use base ('Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf');

## Hive common configuration (every hive pipeline needs this, i.e for using "INPUT_PLUS()")
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

#defines default values for some of the parameters
sub default_options {
  my ($self) = @_;
  return {
    %{$self->SUPER::default_options},
P
    skip_blast_update          => 0, # By default, each blast DB will be checksummed to test for an eventual update (safest). 1 avoids this step saving time.
    min_blast_identity         => 100,  #by default, do not annotate proteins that do not have 100 % identity
    current_query_length_ratio => 1,  #by default, do not annotate the blasted proteins that do not have equal length (ratio 1) than their query
    write                      => 0, #do not write the xrefs to the db unless explicitely being told so.
        
    # the following parameters are part of DBFactory    
    species      => [], 
    taxons       => [],
    division     => [],
    antispecies  => [],
    antitaxons   => [],
    run_all      => 0,
    meta_filters => {},
        
    db_type => 'core',
    unique_rows => 1,
  };
}

#Defines which parameters are required from the user command's line
sub pipeline_wide_parameters {
  my ($self) = @_;
  return {
     %{$self->SUPER::pipeline_wide_parameters},

    'inputfile'             => $self->o('input_file'),
    'blast_db_directory'    => $self->o('blast_db_dir'),    
    'phi_release'           => $self->o('phi_release'),
    '_division'             => $self->o('division'),
  };
}

=head2 pipeline_analyses

    Description : Implements pipeline_analyses() interface method of Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that defines the structure of the pipeline: analyses, jobs, rules, etc.

=cut

sub pipeline_analyses {
  my ($self) = @_;
  
  return [
    {
      -logic_name => 'inputfile',
      -module     => 'Bio::EnsEMBL::EGPipeline::PHIBase::PHIFileReader',
      -input_ids  => [{
                       #seeding the pipeline from user provided value
                       'inputfile' => $self->o('input_file'),
                     }],
      -rc_name    => '4Gb_job',
      -parameters => {
                       delimiter => ',',
                       column_names => 1,
                       output_ids => '#output_ids#',
                       inputfile => '#inputfile#',
                     },
      -flow_into  => {
                       '2->A' => { 'find_subtaxons' => INPUT_PLUS()},
                       'A->1' => ['dbFactory_clean_xrefs'],
                     },
    },
    {
      -logic_name         => 'find_subtaxons',
      -module             => 'Bio::EnsEMBL::EGPipeline::PHIBase::SubtaxonsFinder',
      -analysis_capacity  => 50,
      -batch_size         => 200,
      -max_retry_count    => 2,
      -rc_name            => '4Gb_job',
      -flow_into          => {
                               2 => { 'find_translation' => INPUT_PLUS() },
                             },
    },
    {
      -logic_name         => 'find_translation',
      -module             => 'Bio::EnsEMBL::EGPipeline::PHIBase::TranslationFinder',
      -parameters         => {
                               skip_blast_update => $self->o('skip_blast_update'),
                             },
      -analysis_capacity  => 50,
      -batch_size         => 100,
      -max_retry_count    => 2,
      -flow_into          => {
          3 => WHEN('(#_evidence# eq "SEQUENCE_MATCH") && (#skip_blast_update# == 0)' => { 'blast_p' => INPUT_PLUS()},),
          4 => WHEN ('(#_evidence# eq "SEQUENCE_MATCH") && (#skip_blast_update# == 1)' => { 'blast_p_skip_update' => INPUT_PLUS()},),
          1 => WHEN ( '#_evidence# eq "DIRECT"' => {'phi_accumulator' => INPUT_PLUS() },),
      },
    },
    {
      -logic_name    => 'blast_p',
      -module        => 'Bio::EnsEMBL::EGPipeline::PHIBase::ProteinBlaster',
      -parameters    => {
                          skip_blast_update => $self->o('skip_blast_update'),
                          min_blast_identity => $self->o('min_blast_identity'),
                          current_query_length_ratio => $self->o('current_query_length_ratio'),
                        },
      -rc_name            => '32Gb_job',
      -analysis_capacity  => 100,
      -batch_size         => 200,
      -max_retry_count    => 2,
      -flow_into          => {
                               3 => ['phi_accumulator' ],
                             },
    },
    {
      -logic_name    => 'blast_p_skip_update',
      -module        => 'Bio::EnsEMBL::EGPipeline::PHIBase::ProteinBlaster',
      -parameters    => {
                          skip_blast_update => $self->o('skip_blast_update'),
                          min_blast_identity => $self->o('min_blast_identity'),
                          current_query_length_ratio => $self->o('current_query_length_ratio'),
                          },
      -rc_name       => 'default',
      -analysis_capacity  => 100,
      -batch_size         => 200,
      -max_retry_count    => 2,
      -flow_into          => {
                               3 => ['phi_accumulator' ],
                             },
    },
    {
      -logic_name    => 'phi_accumulator',
      -module        => 'Bio::EnsEMBL::EGPipeline::PHIBase::PHIAccumulatorWriter',
      -analysis_capacity  => 100,
      -batch_size         => 250,
      -max_retry_count    => 2,
      -flow_into => {                    
                      1 => ['?accu_name=phibase_xrefs&accu_address={_species}{translation_db_id}{phi_entry}[]&accu_input_variable=associated_xrefs'],
                    },
    },
    {
      -logic_name         => 'clean_xrefs',
      -module             => 'Bio::EnsEMBL::EGPipeline::PHIBase::XrefsCleaner',
      -parameters         => {
                               append      => [qw(-N)],
                               db_type     => $self->o('db_type'),
                             },
      -analysis_capacity  => 100,
      -max_retry_count    => 1,
      -batch_size         => 1,
      -rc_name            => 'default',
    },            
    {
      -logic_name    => 'dbload_xrefs',
      -module        => 'Bio::EnsEMBL::EGPipeline::PHIBase::XrefsDbLoader',
      -parameters    => {
                          write => $self->o('write'),
                        },
      -flow_into     => {
                          1 => { 'dbFactory_summary' => INPUT_PLUS() },
                        },
    },
    {
      -logic_name      => 'dbFactory_clean_xrefs',
      -module          => 'Bio::EnsEMBL::EGPipeline::PHIBase::DbFactory',
      -max_retry_count => 0,
      -parameters      => {
                            species      => $self->o('species'),
                            taxons       => $self->o('taxons'),
                            division     => $self->o('division'),
                            antispecies  => $self->o('antispecies'),
                            antitaxons   => $self->o('antitaxons'),
                            run_all      => $self->o('run_all'),
                            meta_filters => $self->o('meta_filters'),
                            db_type      => $self->o('db_type'),
                            input_query  => 'DELETE x.*,ox.*,ax.*,ag.*,oox.* 
                                               from external_db e
                                               join xref x using (external_db_id)
                                               join object_xref ox using (xref_id)
                                               join associated_xref ax using (object_xref_id)
                                               join associated_group ag using (associated_group_id)
                                               join ontology_xref oox using (object_xref_id)
                                              WHERE e.db_name="PHI";

                                             DELETE from gene_attrib
                                              WHERE attrib_type_id=358;

                                             DELETE from gene_attrib
                                              WHERE attrib_type_id=317 and value="PHI";',
                          },
      -rc_name         => 'default',
      -flow_into       => {
                            '2->A' => ['clean_xrefs'],
                            'A->1' => ['dbload_xrefs'],
                          },
    },
    {
      -logic_name      => 'dbFactory_summary',
      -module          => 'Bio::EnsEMBL::EGPipeline::PHIBase::DbFactory',
      -max_retry_count => 0,
      -parameters      => {
                            species      => $self->o('species'),
                            taxons       => $self->o('taxons'),
                            division     => $self->o('division'),
                            antispecies  => $self->o('antispecies'),
                            antitaxons   => $self->o('antitaxons'),
                            run_all      => $self->o('run_all'),
                            meta_filters => $self->o('meta_filters'),
                            db_type      => $self->o('db_type'),
                            input_query  => undef,
                          },
      -rc_name         => 'default',
      -flow_into       => {
                            '2->A' => ['phi_summary_sql'],
                            'A->1' => ['build_phi_summary']
                          },
    },
    {
      -logic_name    => 'phi_summary_sql',
      -module        => 'Bio::EnsEMBL::EGPipeline::PHIBase::PHISummarySql',
      -parameters    => {
                          append   => [qw(-N)],
                          db_type  => $self->o('db_type'),
                        },
      -analysis_capacity  => 100,
      -batch_size         => 1,
      -max_retry_count    => 1,
    },
    {
      -logic_name    => 'build_phi_summary',
      -module        => 'Bio::EnsEMBL::EGPipeline::PHIBase::PHISummaryBuilder',
      -parameters    => {
                          append   => [qw(-N)],
                          db_type  => $self->o('db_type'),
                        },
    },
  ];
}

sub hive_meta_table {
  my ($self) = @_;
  return {
    %{$self->SUPER::hive_meta_table},
    'hive_use_param_stack'  => 1, #Jobs see parameters from their ascendants without needing INPUT_PLUS()
  };
}


sub resource_classes {
  my ($self) = @_;
  my $reg_requirement = '--reg_conf ' . $self->o('registry'); #pass registry on to the workers without needing to specify it with the beekeeper

  return {
    %{$self->SUPER::resource_classes},  # inherit 'default' from the parent class
    'default'   => {'LSF' => ['-C0 -q ' . $self->o('queue_name') . ' -M 100   -R"select[mem>100]   rusage[mem=100]"',   $reg_requirement], 'LOCAL' => ['', $reg_requirement] },
    '4Gb_job'   => {'LSF' => ['-C0 -q ' . $self->o('queue_name') . ' -M 4000  -R"select[mem>4000]  rusage[mem=4000]"',  $reg_requirement], 'LOCAL' => ['', $reg_requirement] },
    '8Gb_job'   => {'LSF' => ['-C0 -q ' . $self->o('queue_name') . ' -M 8000  -R"select[mem>8000]  rusage[mem=8000]"',  $reg_requirement], 'LOCAL' => ['', $reg_requirement] },
    '10Gb_job'  => {'LSF' => ['-C0 -q ' . $self->o('queue_name') . ' -M 10000 -R"select[mem>10000] rusage[mem=10000]"', $reg_requirement], 'LOCAL' => ['', $reg_requirement] },
    '16Gb_job'  => {'LSF' => ['-C0 -q ' . $self->o('queue_name') . ' -M 16000 -R"select[mem>16000] rusage[mem=16000]"', $reg_requirement], 'LOCAL' => ['', $reg_requirement] },
    '32Gb_job'  => {'LSF' => ['-C0 -q ' . $self->o('queue_name') . ' -M 32000 -R"select[mem>32000] rusage[mem=32000]"', $reg_requirement], 'LOCAL' => ['', $reg_requirement] },
  };
}

1;
