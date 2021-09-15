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

EGGeneric_conf

=head1 DESCRIPTION

EG specific extensions to the generic hive config.
Serves as a single place to configure EG pipelines.

=cut

=head2 default_options

Description: Interface method that should return a hash of
             option_name->default_option_value pairs.
             
=cut

package perl::Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::Version 2.2;
use perl::Bio::EnsEMBL::EGPipeline::PrivateConfDetails;
use base ('Bio::EnsEMBL::Hive::PipeConfig::EnsemblGeneric_conf');

sub default_options {
  my ($self) = @_;
  return {
    # Inherit options from the base class.
    # Useful ones are:
    #  ensembl_cvs_root_dir (defaults to env. variable ENSEMBL_CVS_ROOT_DIR)
    #  ensembl_release (retrieved from the ApiVersion module)
    #  dbowner (defaults to env. variable USER)
    # The following db-related variables exist in the base class,
    # but will almost certainly need to be overwritten, either by
    # specific *_conf.pm files or on the command line:
    #  host, port, user, password, pipeline_name.
    # These variables are the default parameters used to create the pipeline_db.
    %{$self->SUPER::default_options},
    
    # Generic EG-related options.
    email => $self->o('ENV', 'USER').'@ebi.ac.uk',
    
    # Don't fall over if someone has the temerity to use 'pass' instead of 'password'
    pass => $self->o('password'),
    password => $self->o('pass'),
    
    # Allow a bit of flexibility in the naming of the registry parameter
    reg_conf => $self->o('registry'),
    registry => $self->o('reg_conf'),

    # LinuxBrew home path
    linuxbrew_home  => $self->o('ENV', 'LINUXBREW_HOME'),

    # default LSF queueu name
    queue_name =>  $self->private_conf('ENSEMBL_QUEUE_NAME'),
    datamove_queue_name =>  $self->private_conf('ENSEMBL_DATAMOVE_QUEUE_NAME'),

    # pipeline tag
    pipeline_tag => '',
  }
}

=head2 resource_classes

Description: Interface method that should return a hash of
             resource_description_id->resource_description_hash.
             
=cut

sub resource_classes {
  my ($self) = @_;
  return {
    'datamove'          => {'LSF' => '-q ' . $self->o('datamove_queue_name')},
    'default'           => {'LSF' => '-q standard ' . $self->o('queue_name') . ' -M  4000 -R "rusage[mem=4000]"'},
    'normal'            => {'LSF' => '-q standard ' . $self->o('queue_name') . ' -M  4000 -R "rusage[mem=4000]"'},
    '2Gb_mem'           => {'LSF' => '-q standard ' . $self->o('queue_name') . ' -M  2000 -R "rusage[mem=2000]"'},
    '4Gb_mem'           => {'LSF' => '-q standard ' . $self->o('queue_name') . ' -M  4000 -R "rusage[mem=4000]"'},
    '8Gb_mem'           => {'LSF' => '-q standard ' . $self->o('queue_name') . ' -M  8000 -R "rusage[mem=8000]"'},
    '12Gb_mem'          => {'LSF' => '-q ' . $self->o('queue_name') . ' -M 12000 -R "rusage[mem=12000]"'},
    '16Gb_mem'          => {'LSF' => '-q ' . $self->o('queue_name') . ' -M 16000 -R "rusage[mem=16000]"'},
    '24Gb_mem'          => {'LSF' => '-q ' . $self->o('queue_name') . ' -M 24000 -R "rusage[mem=24000]"'},
    '32Gb_mem'          => {'LSF' => '-q ' . $self->o('queue_name') . ' -M 32000 -R "rusage[mem=32000]"'},
    '2Gb_mem_4Gb_tmp'   => {'LSF' => '-q ' . $self->o('queue_name') . ' -M  2000 -R "rusage[mem=2000,tmp=4000]"'},
    '4Gb_mem_4Gb_tmp'   => {'LSF' => '-q ' . $self->o('queue_name') . ' -M  4000 -R "rusage[mem=4000,tmp=4000]"'},
    '8Gb_mem_4Gb_tmp'   => {'LSF' => '-q ' . $self->o('queue_name') . ' -M  8000 -R "rusage[mem=8000,tmp=4000]"'},
    '12Gb_mem_4Gb_tmp'  => {'LSF' => '-q ' . $self->o('queue_name') . ' -M 12000 -R "rusage[mem=12000,tmp=4000]"'},
    '16Gb_mem_4Gb_tmp'  => {'LSF' => '-q ' . $self->o('queue_name') . ' -M 16000 -R "rusage[mem=16000,tmp=4000]"'},
    '16Gb_mem_16Gb_tmp' => {'LSF' => '-q ' . $self->o('queue_name') . ' -M 16000 -R "rusage[mem=16000,tmp=16000]"'},
    '24Gb_mem_4Gb_tmp'  => {'LSF' => '-q ' . $self->o('queue_name') . ' -M 24000 -R "rusage[mem=24000,tmp=4000]"'},
    '32Gb_mem_4Gb_tmp'  => {'LSF' => '-q ' . $self->o('queue_name') . ' -M 32000 -R "rusage[mem=32000,tmp=4000]"'},
  }
}

=head2 check_exe_in_cellar

Description: Interface method that should return path to executable in Cellar (LinuxBrew home relative).

=cut

# see ensembl Bio::EnsEMBL::Compara::PipeConfig::ComparaGeneric_conf
sub check_exe_in_cellar {
  my ($self, $exe_path) = @_;
  $exe_path = "Cellar/$exe_path";
  push @{$self->{'_all_exe_paths'}}, $exe_path;
  return $self->o('linuxbrew_home').'/'.$exe_path;
}

=head2 tagged_pipeline_name

Description: Interface method that should return pipeline name followed by the 'pipeline_tag' value

=cut

sub tagged_pipeline_name {
  my ($self, @tags) = @_;
  return join('_', grep { defined $_ } @tags) . $self->o('pipeline_tag');
}

sub private_conf {
  my ($self, $conf_key) = @_;
  return Bio::EnsEMBL::EGPipeline::PrivateConfDetails::private_conf($conf_key);
}

1;
