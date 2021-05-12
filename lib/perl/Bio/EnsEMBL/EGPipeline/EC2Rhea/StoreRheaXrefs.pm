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

package Bio::EnsEMBL::EGPipeline::EC2Rhea::StoreRheaXrefs;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');
use Data::Dumper;

use Bio::EnsEMBL::OntologyXref;

sub run {
  my ($self) = @_;

  my $aa   = $self->core_dba->get_adaptor('Analysis');
  my $dbea = $self->core_dba->get_adaptor('DBEntry');
  my $dba = $self->get_DBAdaptor("core");

  my $species_id = $self->core_dba->species_id;
  my $external_db = "Rhea";

  print $species_id . "\n";

  # Get hash of EC numbers to array of translation_ids using SQL

  my $ec_to_translation_ids = $self->fetch_translations_with_xrefs($dba, "IntEnz");

  # Get hash of EC numbers to array of Rhea ids from the file

  my $ec_to_master_rhea = $self->parse_ec2rhea;

  # Delete all existing xrefs for Rhea

  my $external_dbs = [$external_db];

  $self->remove_xrefs($dba, $external_dbs);

  # Get list of all Rhea Xref to store that have EC numbers on your database
  # and construct an xref for each of them, return array of Xref with list of
  # translation_ids

  my $analysis = $aa->fetch_by_logic_name('ec2rhea');

  $self->add_xrefs($dba, $analysis, $external_db, $ec_to_translation_ids, $ec_to_master_rhea);

}

sub add_xrefs {
  my ($self, $dba, $analysis, $external_db, $ec_to_translation_ids, $ec_to_master_rhea) = @_;

  my $ta   = $dba->get_adaptor('Translation');
  my $dbea = $dba->get_adaptor('DBEntry');

  foreach my $ec_string ( keys %$ec_to_translation_ids ) {
    if ($ec_to_master_rhea->{$ec_string}) {
      my @rhea_ids = @ { $ec_to_master_rhea->{$ec_string} };

      foreach my $rhea_id ( @rhea_ids ) {
        my $xref = $self->add_xref($rhea_id, $analysis, $external_db);

        my @list_of_translations = @{ $ec_to_translation_ids->{$ec_string} };

        #For each translation do a STORE
        foreach my $dbID (@list_of_translations) {
          $dbea->store($xref, $dbID, 'Translation', undef);
        }
      }
    }
  }
}

sub add_xref {
  my ($self, $rhea_id, $analysis, $external_db) = @_;

  my $xref = Bio::EnsEMBL::DBEntry->new(
    -PRIMARY_ID => $rhea_id,
    -DISPLAY_ID => $rhea_id,
    -DBNAME     => $external_db,
    -INFO_TYPE  => 'PROJECTION',
  );
  $xref->analysis($analysis);

  return $xref;
}

sub fetch_translations_with_xrefs {
  my ($self, $dba, $external_db_name) = @_;

  my $sql = q/
    SELECT x.display_label, tn.translation_id FROM
      translation tn INNER JOIN
      transcript tt USING (transcript_id) INNER JOIN
      seq_region sr USING (seq_region_id) INNER JOIN
      coord_system cs USING (coord_system_id) INNER JOIN
      object_xref ox ON (ox.ensembl_id = tn.translation_id) INNER JOIN
      xref x USING (xref_id) INNER JOIN
      external_db edb USING (external_db_id) LEFT OUTER JOIN
      identity_xref ix USING (object_xref_id) LEFT OUTER JOIN
      ontology_xref ontix USING (object_xref_id)
    WHERE
      ox.ensembl_object_type = 'Translation' AND
      cs.species_id = ? AND
      edb.db_name = ?
  /;

  my $sth = $dba->dbc->prepare($sql);
  $sth->execute($dba->species_id, $external_db_name);

  my $results = $sth->fetchall_arrayref();

  my %ec_to_translation_ids;

  foreach my $line ( @{$results} ) {
    push( @{ $ec_to_translation_ids { $line->[0] } }, $line->[1]) unless $line->[0] =~ /-/;
  }

  return \%ec_to_translation_ids;

}

sub parse_ec2rhea {
  my ($self) = @_;

  my $ec2rhea_file = $self->param_required('ec2rhea_file');

  open (my $fh, '<', $ec2rhea_file) or die "Failed to open $ec2rhea_file: $!\n";

  my %ec_to_rhea;
  while (my $line = <$fh>) {
    next if ($line =~ /^RHEA_ID/);

    my ($rhea_id, $direction, $master_id, $ec_id) = $line
      =~ /^(\d+)\s+(\w+)\s+(\d+)\s+(\S+)/;

    push( @{ $ec_to_rhea {$ec_id} }, $master_id);

  }

  close $fh;

  return \%ec_to_rhea;
}

sub remove_xrefs {
  my ($self, $dba, $external_dbs) = @_;

  my $db_names = "'" . join("','", @$external_dbs) . "'";

  my $sql = q/
    DELETE ox.*, ix.*, ontix.* FROM
      translation tn INNER JOIN
      transcript tt USING (transcript_id) INNER JOIN
      seq_region sr USING (seq_region_id) INNER JOIN
      coord_system cs USING (coord_system_id) INNER JOIN
      object_xref ox ON (ox.ensembl_id = tn.translation_id) INNER JOIN
      xref x USING (xref_id) INNER JOIN
      external_db edb USING (external_db_id) LEFT OUTER JOIN
      identity_xref ix USING (object_xref_id) LEFT OUTER JOIN
      ontology_xref ontix USING (object_xref_id)
    WHERE
      ox.ensembl_object_type = 'Translation' AND
      cs.species_id = ? AND
  /;
  $sql .= "edb.db_name IN ($db_names)";

  my $sth = $dba->dbc->db_handle->prepare($sql);
  $sth->execute($dba->species_id);

}

1;
