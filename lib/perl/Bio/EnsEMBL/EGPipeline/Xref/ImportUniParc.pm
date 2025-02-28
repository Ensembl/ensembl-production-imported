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

Bio::EnsEMBL::EGPipeline::Xref::ImportUniParc

=head1 DESCRIPTION

Add UniParc checksums to local database.

=head1 Author

James Allen

=cut

use strict;
use warnings;
use feature 'say';

package Bio::EnsEMBL::EGPipeline::Xref::ImportUniParc;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use File::Copy qw(move);
use File::Path qw(make_path);
use File::Spec::Functions qw(catdir);
use File::Temp qw(tempdir);
use Time::Local;

sub param_defaults {
  my ($self) = @_;
  return {
    'ftp_file' => $self->private_conf('ENSEMBL_UNIPARC_UPI_LIST'),
    'tmp_dir'  => $self->private_conf('ENSEMBL_TMPDIR') || '/tmp',
  };
}

# uniparc_dbm_cache_dir => $self->o('uniparc_dbm_cache_dir'),
# uniparc_dbm_cache_name => $self->o('uniparc_dbm_cache_name'),
# dbm_create_script => catdir($self->o('ensembl_production_imported_scripts_dir'), 'uniparc_index', 'create_uniparc_dbm.py'),

sub run {
  my ($self) = @_;
  my $ftp_file  = $self->param_required('ftp_file');
  my $dbm_dir   = $self->param_required('uniparc_dbm_cache_dir');
  my $dbm_name  = $self->param_required('uniparc_dbm_cache_name');
  my $create_py = $self->param_required('dbm_create_script');

  if (! -e $dbm_dir) {
    make_path($dbm_dir) or $self->throw("Failed to create directory '$dbm_dir': $!");
  }

  my $lock_file = catdir($dbm_dir, "lock");
  $self->lock($lock_file);

  my $dbm_file = catdir($dbm_dir, $dbm_name);
  if ($self->import_required($ftp_file, $dbm_file)) {
    $self->import_uniparc($ftp_file, $create_py, $dbm_dir, $dbm_name);
  }

  $self->unlock($lock_file);
}

sub lock {
  my ($self, $lock_file) = @_;
  # locking on NFS!!!
}

sub unlock {
  my ($self, $lock_file) = @_;
  # locking on NFS!!!
}

sub import_required {
  my ($self, $ftp_file, $dbm_file) = @_;
  my $required = 1;

  if (! -f $dbm_file ) {
    return $required;
  }

  my $ftp_file_timestamp = (stat $ftp_file)[9];
  my $dbm_file_timestamp = (stat $dbm_file)[9];

  if ($dbm_file_timestamp > $ftp_file_timestamp) {
    $required = 0;
  }

  return $required;
}

sub import_uniparc {
  my ($self, $ftp_file, $create_py, $dbm_dir, $dbm_name) = @_;

  my $cat= "cat";
  if ($ftp_file =~ m/\.gz$/) {
    $cat = "zcat";
  }

  my $shm_dir = tempdir( DIR => catdir("", "dev", "shm"), CLEANUP => 1 ); # /dev/shm -- os specific, but fast
  my $shm_file = catdir($shm_dir, $dbm_name);
  my $log_file = catdir($dbm_dir, $dbm_name . ".create.log");
  my $log_file_new = catdir($dbm_dir, $dbm_name . ".create.new.log");

  my $cmd = "";
  $cmd .= "$cat $ftp_file | ";
  $cmd .= "  python $create_py --dbfile $shm_file 2> $log_file_new";
  system($cmd) == 0 or $self->throw("Failed to run '$cmd': $!");

  move($shm_file, $dbm_dir) or $self->throw("Failed to move '$shm_file' to '$dbm_dir': $!");
  move($log_file_new, $log_file) or $self->throw("Failed to move '$log_file_new' to '$log_file': $!");
}

1;
