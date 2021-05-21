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
use File::Copy qw(copy);
use Time::Local;

sub param_defaults {
  my ($self) = @_;
  return {
    'ftp_file' => $self->private_conf('ENSEMBL_UNIPARC_UPI_LIST'),
    'tmp_dir'  => $self->private_conf('ENSEMBL_TMPDIR') || '/tmp',
  };
}

sub run {
  my ($self) = @_;
  my $ftp_file    = $self->param_required('ftp_file');
  my $tmp_dir     = $self->param_required('tmp_dir');
  my $uniparc_db  = $self->param_required('uniparc_db');
  my $uniparc_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%$uniparc_db);

  if ($self->import_required($ftp_file, $uniparc_dba)) {
    $self->import_uniparc($ftp_file, $tmp_dir, $uniparc_dba);
  }
}

sub import_required {
  my ($self, $ftp_file, $uniparc_dba) = @_;
  my $required = 1;

  my $file_timestamp = (stat $ftp_file)[9];

  my $uniparc_dbh = $uniparc_dba->dbc->db_handle();

  my $sql = "SELECT update_time FROM last_update ORDER BY update_time DESC";
  my $sth = $uniparc_dbh->prepare($sql);
  $sth->execute();

  my $results = $sth->fetchall_arrayref();
  if (scalar(@$results)) {
    my $update_time = $$results[0][0];
    my ($year, $month, $day, $h, $m, $s) = split(/\D/, $update_time);
    my $update_time_epoch = timelocal($s, $m, $h, $day, $month-1, $year-1900);

    if ($update_time_epoch > $file_timestamp) {
      $required = 0;
    }
  }

  return $required;
}

sub import_uniparc {
  my ($self, $ftp_file, $tmp_dir, $uniparc_dba) = @_;

  copy $ftp_file, "$tmp_dir/protein.txt";

  my $uniparc_dbh = $uniparc_dba->dbc->db_handle();

  my @preload_sql = (
    "DROP INDEX md5_idx ON protein;",
    "TRUNCATE TABLE protein;",
  );

  foreach my $load_sql (@preload_sql) {
    $uniparc_dbh->do($load_sql) or $self->throw("Failed to execute: $load_sql");
  }

  my $cmd = $self->mysqlimport_command_line($uniparc_dba->dbc);
  $cmd .= " --fields_terminated_by ' ' $tmp_dir/protein.txt";
  system($cmd) == 0 or $self->throw("Failed to run ".$cmd);

  my @postload_sql = (
    "CREATE INDEX md5_idx ON protein(md5);",
    "TRUNCATE TABLE last_update;",
    "INSERT INTO last_update (update_time) VALUES (NOW());",
  );

  foreach my $load_sql (@postload_sql) {
    $uniparc_dbh->do($load_sql) or $self->throw("Failed to execute: $load_sql");
  }
}

1;
