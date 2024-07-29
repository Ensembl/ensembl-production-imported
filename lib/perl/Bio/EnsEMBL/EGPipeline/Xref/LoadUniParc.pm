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

Bio::EnsEMBL::EGPipeline::Xref::LoadUniParc

=head1 DESCRIPTION

Add UniParc xrefs to a core database, based on checksums.

=head1 Author

Dan Staines and James Allen

=cut

package Bio::EnsEMBL::EGPipeline::Xref::LoadUniParc;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::Xref::LoadXref');

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Digest::MD5;

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    'logic_name'  => 'xrefchecksum',
    'external_db' => 'UniParc',
  };
}

sub run {
  my ($self) = @_;
  my $db_type     = $self->param_required('db_type');
  my $logic_name  = $self->param_required('logic_name');
  my $external_db = $self->param_required('external_db');

  my $dba = $self->get_DBAdaptor($db_type);
  my $aa  = $dba->get_adaptor('Analysis');

  my $analysis = $aa->fetch_by_logic_name($logic_name);

  $self->external_db_reset($dba, $external_db);

  $self->add_xrefs($dba, $analysis, $external_db);

  $self->external_db_update($dba, $external_db);

  $dba->dbc and $dba->dbc->disconnect_if_idle();
}

sub add_xrefs {
  my ($self, $dba, $analysis, $external_db) = @_;
  my $uniparc_db = $self->param_required('uniparc_db');

  my $uniparc_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%$uniparc_db);

  my $ta   = $dba->get_adaptor('Translation');
  my $dbea = $dba->get_adaptor('DBEntry');

  my $translations = $ta->fetch_all();
  foreach my $translation (@$translations) {
    my $upi = $self->search_for_upi($uniparc_dba, $translation);
    if ($upi) {
         my $xref = $self->add_xref($upi, $analysis, $external_db);
          $dbea->store($xref, $translation->dbID(), 'Translation');
    }
  }
}

sub add_xref {
  my ($self, $upi, $analysis, $external_db) = @_;

  my $xref = Bio::EnsEMBL::DBEntry->new(
    -PRIMARY_ID => $upi,
    -DISPLAY_ID => $upi,
    -DBNAME     => $external_db,
    -INFO_TYPE  => 'CHECKSUM',
  );
  $xref->analysis($analysis);

  return $xref;
}

sub search_for_upi {
  my ($self, $uniparc_dba, $translation) = @_;

  my $checksum = $self->md5_checksum($translation->seq);

  my $sql = 'SELECT upi FROM protein WHERE md5 = ?';
  my $sth = $uniparc_dba->dbc->db_handle->prepare($sql);
  eval { $sth->execute($checksum) };
  if ($sth->err) {
    $self->warning("Unable to search for UPI(s) of protein with MD5 $checksum, DB error: $sth->err : $sth->errstr");
    return;
  }

  my $upi;
  my $results = $sth->fetchall_arrayref();
  if (scalar(@$results)) {
    if (scalar(@$results) == 1) {
      $upi = $$results[0][0];
    } else {
      $self->warning("Multiple UPIs found for ".$translation->stable_id);
    }
  } else {
    $self->warning("No UPI found for ".$translation->stable_id);
  }

  return $upi;
}

sub md5_checksum {
  my ($self, $sequence) = @_;

  my $digest = Digest::MD5->new();
  if ($sequence =~ /^X[^X]/) {
    $sequence =~ s/^X//;
  }
  $digest->add($sequence);

  return uc($digest->hexdigest());
}

1;
