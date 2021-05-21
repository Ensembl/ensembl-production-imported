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

package Bio::EnsEMBL::EGPipeline::SequenceAlignment::Exonerate::StartServer;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub param_defaults {
  my $self = shift @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'tmp_dir'         => $self->private_conf('ENSEMBL_TMPDIR') || '/tmp',
    'log_file'        => $self->o('tmp_dir') . '/exonerate-server.out',
    'max_connections' => 10,
  };
}

sub run {
  my $self = shift @_;
  
  my $server_file = $self->param_required('server_file');
  
  my $server = $self->start_server;
  open my $fh, '>', $server_file or $self->throw("Failed to create $server_file");
  print $fh $server;
  close $fh;
  
  while (-e $server_file) {
    sleep 60;
  }
  
  $self->stop_server;
}

sub start_server {
  my $self = shift @_;
  
  my $server_exe = $self->param_required('server_exe');
  my $index_file = $self->param_required('index_file');
  my $max_connections = $self->param('max_connections');
  
  # Start at default port; if something is already running, do a limited
  # number of increments to see if we can find an available port.
  my $port = 12886;
  my $used = 1;
  while ($used && $port < 12896) {
    $used = `netstat -an | grep :$port`;
    $port++ if $used;
  }
  if ($used) {
    $self->throw("Failed to find an available port for exonerate-server");
  }
  my $log_file = $self->param_required('log_file').".$port";
  $self->param('log_file', $log_file);
  
  my $command = "$server_exe $index_file --maxconnections $max_connections --port $port &> $log_file";
  
  my $pid;
  {
    if ($pid = fork) {
      last;
    } elsif (defined $pid) {
      exec("exec $command") == 0 or $self->throw("Failed to run $command: $!");
    }
  }
  
  my ($server_starting, $cycles) = (1, 0);
  while ($server_starting) {
    if ($cycles < 50) {
      sleep 30;
      $cycles++;
      my $started_message = `tail -1 $log_file`;
      if ($started_message =~ /Message: listening on port/) {
        $server_starting = 0;
      }
    } else {
      $self->throw("Failed to start server; see log: $log_file");
    }
  }
  
  $self->param('server_pid', $pid);
  
  return $self->worker->meadow_host.":$port";
}

sub stop_server {
  my $self = shift @_;
  
  my $pid = $self->param('server_pid');
  $self->warning("Killing server process $pid");
  kill('KILL', $pid) or $self->throw("Failed to kill server process $pid: $!");
}

1;
