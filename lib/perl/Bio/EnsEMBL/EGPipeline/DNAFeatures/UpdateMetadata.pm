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

package Bio::EnsEMBL::EGPipeline::DNAFeatures::UpdateMetadata;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub run {
  my ($self) = @_;
  my $species = $self->param_required('species');
  
  my $dba = $self->core_dba();
  my $dbh = $dba->dbc->db_handle();
  
  # Add meta data
  my $meta_sql =
    "INSERT IGNORE INTO meta (species_id, meta_key, meta_value) ".
    "SELECT " . $dba->species_id . ", 'repeat.analysis', logic_name FROM ".
    "repeat_feature INNER JOIN analysis USING (analysis_id) ".
    "GROUP BY logic_name;";
  $dbh->do($meta_sql);
  
  # Update repeat consensus entries
  my $rc_sql =
    "UPDATE repeat_consensus SET repeat_type = ? ".
    "WHERE repeat_class LIKE ?;";
  my $sth = $dbh->prepare($rc_sql);
  
  my %mappings = (
    'Low_Comp%'  => 'Low complexity regions',
		'LINE%'      => 'Type I Transposons/LINE',
		'SINE%'      => 'Type I Transposons/SINE',
		'DNA%'       => 'Type II Transposons',
		'LTR%'       => 'LTRs',
		'Other%'     => 'Other repeats',
		'Satelli%'   => 'Satellite repeats',
		'Simple%'    => 'Simple repeats',
		'Other%'     => 'Other repeats',
		'Tandem%'    => 'Tandem repeats',
		'TRF%'       => 'Tandem repeats',
		'Waterman'   => 'Waterman',
		'Recon'      => 'Recon',
		'MaskRegion' => 'Mask region',
		'dust%'      => 'Dust',
		'Unknown%'   => 'Unknown',
		'%RNA'       => 'RNA repeats',
  );
	foreach (keys %mappings) {
    $sth->execute($mappings{$_}, $_);
	}
  
  my $rc_unknown_sql =
    "UPDATE repeat_consensus SET repeat_type = 'Unknown' ".
    "WHERE repeat_type = '' OR repeat_type IS NULL;";
  $dbh->do($rc_unknown_sql);
  
}

1;
