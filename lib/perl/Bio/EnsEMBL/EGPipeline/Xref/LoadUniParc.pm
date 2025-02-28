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
use File::Path qw(make_path);
use File::Spec::Functions qw(catdir);
use File::Touch;
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

  my $dbm_file = $self->param_required('uniparc_dbm_cache');
  my $query_py = $self->param_required('dbm_query_script');
  my $work_dir = $self->param_required('upi_query_dir'); 

  # dump translation dbIDs, stable IDs, DM5(upper case), sequence
  my $ta = $dba->get_adaptor('Translation');

  if (! -e $work_dir) {
    make_path($work_dir) or $self->throw("Failed to create directory '$work_dir': $!");
  }

  # get translations
  touch( catdir($work_dir, "get_translations_start.tag") );
  my $translations = []; # [ {dbId, stableID, md5, upis, sequence} ]
  $self->get_translations($ta, $translations);
  touch( catdir($work_dir, "get_translations_end.tag") );

  my $translations_raw_tsv = catdir($work_dir, "translations.raw.tsv");
  $self->dump_translations($translations, $translations_raw_tsv);

  # get unique md5s and query fo UPIs
  my $unique_md5s = {};
  $self->get_uniq_md5s($translations, $unique_md5s);

  my $known_upis = {}; # md5 -> [ upis ]
  $self->query_md5_cache($unique_md5s, $known_upis);

  # dump result
  my $translations_upi_tsv = catdir($work_dir, "translations.upi.tsv");
  $self->dump_translations($translations, $translations_raw_tsv);

  # store xrefs
  touch( catdir($work_dir, "xrefs_store_start.tag") );
  my $translations_store_log = catdir($work_dir, "translations_store.log");
  $self->store_xrefs($dba, $translations, $translations_store_log);
  touch( catdir($work_dir, "xrefs_store_stop.tag") );

  die "";
}


sub store_xrefs {
  my ($self, $dba, $translations, $log) = @_;

  open(my $fh, ">", $log)
    or $self->throw("Failed to open '$log': $!");

  my $dbea = $dba->get_adaptor('DBEntry');

  foreach my $translation (@$translations) {
    my $dbID = $translation->{dbID};
    my $stable_id = $translation->{stable_id} // "";
    my $md5 = $translation->{md5};
    my $upis = $translation->{upis};

    my $upis_cnt = scalar(@$upis);
    if ($upis_cnt > 0) {
      if ($upis_cnt > 1) {
        my $upis_str = join(",", @{$upis});
        my $warnig_str = "Multiple UPIs ($upis_str) found for translation $dbID ($stable_id) with MD5 $md5";
        $self->warning($warning_str);
        print $fh "$warning_str\n";
      }
      my $upi = $upis->[0];
      print $fh "Storing xref $upi for translation $dbID ($stable_id) with MD5 $md5\n";
      #my $xref = $self->add_xref($upis->[0], $analysis, $external_db);
      #$dbea->store($xref, $translation->{dbID}, 'Translation');
    }
  }

  close($fh);
}

sub get_translations {
  my ($self, $ta, $store) = @_;

  my $translations = $ta->fetch_all();
  foreach my $translation (@$translations) {
    my $seq = $translation->seq;
    my $checksum = $self->md5_checksum($seq);
    my $item = {
      dbID => $translation->dbID(),
      stable_id => $translation->stable_id,
      sequence => $seq,
      md5 => $checksum,
      upis => [],
    };

    push @$store, $item;
}

sub dump_translations {
  my ($self, $translations, $file) = @_;

  return if (!file);

  open(my $fh, ">", $file)
    or $self->throw("Failed to open '$file': $!");

  foreach my $item (@$translations) {
    my $upis_str = join(",", @{$item->{upis}});
    my $item_str = join("\t", $item->{dbID}, $item->{stable_id} // "", $item->{md5}, $upis_str, $item->{sequence});
    print $fh "$item_str\n";
  }
  
  close($fh);
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
