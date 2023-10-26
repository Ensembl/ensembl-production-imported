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

package Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::Version 2.2;
use Bio::EnsEMBL::EGPipeline::PrivateConfDetails;
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

    # default LSF queue name
    queue_name =>  $self->private_conf('ENSEMBL_QUEUE_NAME'),
    datamove_queue_name =>  $self->private_conf('ENSEMBL_DATAMOVE_QUEUE_NAME'),

    # pipeline tag
    pipeline_tag => '',
  }
}

=head2 make_resources

Description: Method to generate a resource class for all available meadows (LSF, SLURM), instead of hardcopding the resources.

It should return a hash like { LSF => '', SLURM => '' }

Args: A hash with the following keys:

=over

=item * queue = name of the queue/partition to use

=item * memory = integer in MB

=item * time = time limit in the form 0:00:00 (h:mm:ss)

=item * cpus = number of cores

=item * temp_memory = reserve this amount of temp memory

=item * lsf_param = a string for specific LSF parameters

=item * slurm_param = a string for specific SLURM parameters

=back
             
=cut


sub make_resource {
  my ($self, $conf) = @_;

  my $res = {
    'LSF' => _lsf_resource($conf),
    'SLURM' => _slurm_resource($conf),
  };
  return $res;
}

sub _lsf_resource {
  my ($conf) = @_;

  my $mem = $conf->{memory};
  my $tmem = $conf->{temp_memory};
  my $cpus = $conf->{cpus};
  my $queue = $conf->{queue};
  my $time = $conf->{time};

  if ($time and $time =~ /(\d+):(\d\d):(\d\d)/) {
    $time = "$1:$2";
  }

  my @rusage = ();
  push @rusage, "mem=$mem" if $mem;
  push @rusage, "tmp=$tmem" if $tmem;
  
  my @res_params;
  push @res_params, "-M $mem" if $mem;
  push @res_params, ('-R "rusage[' . join(',', @rusage) . ']"') if @rusage;
  push @res_params, "-q $queue" if $queue;
  push @res_params, "-We $time" if $time;
  push @res_params, "-n $cpus" if $cpus;
  push @res_params, $conf->{lsf_params} if $conf->{lsf_params};
  my $res_string = join(" ", @res_params);
  return $res_string;
}

sub _slurm_resource {
  my ($conf) = @_;

  my $mem = $conf->{memory};
  my $cpus = $conf->{cpus};
  my $queue = $conf->{queue};
  my $time = $conf->{time};

  # Prepare memory string
  my $rmem;
  if ($mem) {
    if ($mem > 1000) {
      $mem = int($mem/1000);
      $rmem = $mem . 'g';
    } else {
      $rmem = $mem . "m";
    }
  }
  
  my @res_params;
  push @res_params, "--mem=$rmem" if $rmem;
  push @res_params, "--time=$time" if $time;
  push @res_params, "-c $cpus" if $cpus;
  push @res_params, "--partition=$queue" if $queue;
  push @res_params, $conf->{slurm_params} if $conf->{slurm_params};
  my $res_string = join(" ", @res_params);
  return $res_string;
}

=head2 resource_classes

Description: Interface method that should return a hash of
             resource_description_id->resource_description_hash.
             
=cut

sub resource_classes {
  my ($self) = @_;

  my $data_queue = $self->o('datamove_queue_name');
  my $queue = $self->o('queue_name');
  my $short = "1:00:00";
  my $long = "24:00:00";

  my %resources = (
    'default'           => $self->make_resource({"queue" => $queue, "memory" => 4_000, "time" => $short}),
    'normal'            => $self->make_resource({"queue" => $queue, "memory" => 4_000, "time" => $long}),
    'datamove'          => $self->make_resource({"queue" => $data_queue, "memory" => 100, "time" => $short}),
    'datamove_4Gb_mem'  => $self->make_resource({"queue" => $data_queue, "memory" => 4_000, "time" => $long}),
    'datamove_32Gb_mem' => $self->make_resource({"queue" => $data_queue, "memory" => 32_000, "time" => $long}),
  );

  my @mems = (2, 4, 8, 12, 16, 32);
  my $tmem = 4;
  my $time = $long;

  for my $mem (@mems) {
    my $name = "${mem}Gb_mem";
    my $tname = "${name}_${tmem}Gb_tmp";

    $resources{$name} = $self->make_resource({"queue" => $queue, "memory" => $mem * 1000, "time" => $time});
    $resources{$tname} = $self->make_resource({"queue" => $queue, "memory" => $mem * 1000, "time" => $time, "temp_memory" => $tmem * 1000});
  }
  $resources{'16Gb_mem_16Gb_tmp'} = $self->make_resource({"queue" => $queue, "memory" => 16_000, "time" => $long, "temp_memory" => 16_000});

  return \%resources;
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
