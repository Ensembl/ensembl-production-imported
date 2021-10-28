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

  perl::Bio::EnsEMBL::EGPipeline::PipeConfig::MapPHI_base_conf

=head1 SYNOPSIS


  init_pipeline.pl perl::Bio::EnsEMBL::EGPipeline::PipeConfig::MapPHI_base_conf  -pipeline_url $EHIVE_URL  -registry $REGISTRY_FILE -blast_db_dir $BLAST_DB_DIRECTORY -input_file $INPUT_ROTHAMSTED_FILE -hive_force_init 1

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

package perl::Bio::EnsEMBL::EGPipeline::PipeConfig::MapPHI_base_conf;


use strict;
use warnings;

## EG common configuration (mostly resource classes)
use base ('perl::Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf');

## Hive common configuration (every hive pipeline needs this, i.e for using "INPUT_PLUS()")
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

#defines default values for some of the parameters
sub default_options {
  my ($self) = @_;
  return {
    %{$self->SUPER::default_options},

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

    'inputfile'             => $self->o('inputfile'),
    'interactions_db_url'   => $self->o('interactions_db_url'),
    'core_db_url'	    => $self->o('core_db_url'), 
    'registry'		    => $self->o('reg_file')
  };
}

=head2 pipeline_analyses

    Description : Implements pipeline_analyses() interface method of Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that defines the structure of the pipeline: analyses, jobs, rules, etc.

=cut

sub pipeline_analyses {
  my ($self) = @_;
  
  return [
    {
      -logic_name => 'input_file',
      -module     => 'ensembl.microbes.runnable.PHIbase_2.FileReader',
      -language   => 'python3',
      -input_ids  => [{
                       #seeding the pipeline from user provided value
                       'inputfile' => $self->o('inputfile'),
		       'registry' => $self->o('reg_file'),
                     }],
      -parameters => {
		       delimiter => ',',
		       column_names => 1,
		       output_ids => '#output_ids#',
		       inputfile => '#inputfile#',
		       registry   => '#registry#',
                      },
      -flow_into    => {
	                2 => { 'meta_ensembl_reader' => INPUT_PLUS() },
		       },
    },
    { 
      -logic_name => 'meta_ensembl_reader',
      -module     => 'ensembl.microbes.runnable.PHIbase_2.MetaEnsemblReader',
      -language   => 'python3',
      -flow_into    => {
                        1 => { 'ensembl_core_reader' => INPUT_PLUS() },
                        -1 => ['failed_entries'],
			},
    },
    {  
       -logic_name => 'ensembl_core_reader',
       -module     => 'ensembl.microbes.runnable.PHIbase_2.EnsemblCoreReader',
       -language   => 'python3',
       -flow_into    => {
                        1 => { 'sequence_finder' => INPUT_PLUS() },
                        -1 => ['failed_entries'],
			},
    },
    { 
       -logic_name => 'sequence_finder',
       -module     => 'ensembl.microbes.runnable.PHIbase_2.SequenceFinder',
       -language   => 'python3',
       -flow_into    => {
                        1  => { 'interactor_data_manager' => INPUT_PLUS() },
			-1 => ['failed_entries'],
                        },
    },
    {
      -logic_name => 'interactor_data_manager',
      -module     => 'ensembl.microbes.runnable.PHIbase_2.InteractorDataManager',
      -language   => 'python3',
      -flow_into    => {
                        1  => ['db_writer'],
			-1 => ['failed_entries'],
                        },
    },
    { 
      -logic_name => 'failed_entries',
      -module     => 'ensembl.microbes.runnable.PHIbase_2.FailedEntries',
      -language   => 'python3',
    },
    { 
      -logic_name => 'db_writer',
      -module     => 'ensembl.microbes.runnable.PHIbase_2.DBwriter',
      -language   => 'python3',
    }
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
  my $reg_requirement = '--reg_conf ' . $self->o('reg_file'); #pass registry on to the workers without needing to specify it with the beekeeper
  return {
    %{$self->SUPER::resource_classes},  # inherit 'default' from the parent class
	     'datamove'          => {'LSF' => '-q ' . $self->o('datamove_queue_name')},
	     'default'           => {'LSF' => '-q ' . $self->o('queue_name') . ' -M  4000 -R "rusage[mem=4000]"'},
	     'normal'            => {'LSF' => '-q ' . $self->o('queue_name') . ' -M  4000 -R "rusage[mem=4000]"'},
	     '2Gb_mem'           => {'LSF' => '-q ' . $self->o('queue_name') . ' -M  2000 -R "rusage[mem=2000]"'},
	     '4Gb_mem'           => {'LSF' => '-q ' . $self->o('queue_name') . ' -M  4000 -R "rusage[mem=4000]"'},
	     '8Gb_mem'           => {'LSF' => '-q ' . $self->o('queue_name') . ' -M  8000 -R "rusage[mem=8000]"'},
	     '12Gb_mem'          => {'LSF' => '-q ' . $self->o('queue_name') . ' -M 12000 -R "rusage[mem=12000]"'},
	     '16Gb_mem'          => {'LSF' => '-q ' . $self->o('queue_name') . ' -M 16000 -R "rusage[mem=16000]"'},
	     '24Gb_mem'          => {'LSF' => '-q ' . $self->o('queue_name') . ' -M 24000 -R "rusage[mem=24000]"'},
	     '32Gb_mem'          => {'LSF' => '-q ' . $self->o('queue_name') . ' -M 32000 -R "rusage[mem=32000]"'},
     };

}

1;

