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
  $self->warning("Getting translations");
  $self->touch( catdir($work_dir, "get_translations_start.tag") );
  my $translations = []; # [ {dbId, stableID, md5, upis, sequence} ]
  $self->get_translations($ta, $translations);
  $self->touch( catdir($work_dir, "get_translations_end.tag") );

  my $translations_raw_tsv = catdir($work_dir, "translations.raw.tsv");
  $self->warning("Dumping raw translations data with no UPIs to $translations_raw_tsv");
  $self->dump_translations($translations, $translations_raw_tsv);

  # query cache and fill UPIs
  $self->warning("Querying and filling UPIs");
  $self->query_cache_and_update_upis($translations, $query_py, $dbm_file, $work_dir);

  # dump result
  my $translations_upi_tsv = catdir($work_dir, "translations.upi.tsv");
  $self->warning("Dumping raw translations data with available UPIs to $translations_upi_tsv");
  $self->dump_translations($translations, $translations_upi_tsv);

  # store xrefs
  $self->touch( catdir($work_dir, "xrefs_store_start.tag") );
  my $translations_store_log = catdir($work_dir, "translations_store.log");
  $self->warning("Storing xrefs. See $translations_store_log");
  $self->store_xrefs($dba, $translations, $analysis, $external_db, $translations_store_log);
  $self->touch( catdir($work_dir, "xrefs_store_stop.tag") );

  $self->warning("Done");
}

sub query_cache_and_update_upis {
  my ($self, $translations, $query_py, $dbm_file, $work_dir) = @_;

  # get unique md5s and query fo UPIs
  $self->warning("Getting unique MD5s");
  my $unique_md5s = [];
  $self->get_uniq_md5s($translations, $unique_md5s);

  # preparing queries
  my $md5_queries = catdir($work_dir, "md5_queries.lst");
  $self->warning("Filling list of MD5 queries to $md5_queries");
  open(my $qfh, ">", $md5_queries)
    or $self->throw("Failed to open '$md5_queries': $!");
  for my $md5 (@$unique_md5s) {
    print $qfh "$md5\n";
  }
  close($qfh);

  # run query script
  my $upi_file = catdir($work_dir, "md5_upi.lst");
  my $log_file = catdir($work_dir, "query.log");

  $self->warning("Fetching UPIs from $dbm_file into $upi_file (log $log_file)");
  my $cmd = "";
  $cmd .= "cat $md5_queries | ";
  $cmd .= "  python $query_py --dbfile $dbm_file > $upi_file 2> $log_file";
  system($cmd) == 0 or $self->throw("Failed to run '$cmd': $!");

  # extract UPIs, assuming same result for the same MD5
  $self->warning("Extracting UPIs from $upi_file");
  my $known_upis = {}; # md5 -> [ upis ]
  open(my $upifh, "<:encoding(utf-8)", $upi_file)
    or $self->throw("Failed to open '$upi_file': $!");
  while(my $line = <$upifh>) {
    chomp $line;
    my ($md5, @upis) = split /\t/, $line;
    $known_upis->{$md5} = [ @upis ];
  }
  close($upifh);

  # update transcript UPIs
  my $known_upis_cnt = scalar(keys(%$known_upis));
  $self->warning("Updating transcripts with $known_upis_cnt known UPIs");
  $self->add_upis($translations, $known_upis);
}

sub get_uniq_md5s {
  my ($self, $translations, $store) = @_;
  push @$store,
    sort { $a cmp $b }
      keys %{{
        map { $_->{md5} => 1 } @$translations
      }};
};

sub add_upis {
  my ($self, $translations, $known) = @_;
  foreach my $translation (@$translations) {
    my $md5 = $translation->{md5};
    my $upis = $translation->{upis};
    if (exists $known->{$md5} && @{$known->{$md5}} ) {
      push @$upis, @{$known->{$md5}};
    }
  }
}

sub store_xrefs {
  my ($self, $dba, $translations, $analysis, $external_db, $log) = @_;

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
        my $warning_str = "Multiple UPIs ($upis_str) found for translation $dbID ($stable_id) with MD5 $md5";
        $self->warning($warning_str);
        print $fh "$warning_str\n";
      }
      my $upi = $upis->[0];
      print $fh "Storing xref $upi for translation $dbID ($stable_id) with MD5 $md5\n";
      my $xref = $self->add_xref($upi, $analysis, $external_db);
      $dbea->store($xref, $dbID, 'Translation');
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
    my $store_size = scalar(@$store);
    if ($store_size % 5_000 == 0) {
      $self->warning("Fetched $store_size translations");
    }
  }
  my $store_size = scalar(@$store);
  $self->warning("Fetched $store_size translations");
}

sub dump_translations {
  my ($self, $translations, $file) = @_;

  return if (!$file);

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

sub md5_checksum {
  my ($self, $sequence) = @_;

  my $digest = Digest::MD5->new();
  if ($sequence =~ /^X[^X]/) {
    $sequence =~ s/^X//;
  }
  $digest->add($sequence);

  return uc($digest->hexdigest());
}

sub touch {
  my ($self, @files) = @_;
  for my $file (@files) {
    system("touch $file");
  }
}


1;
