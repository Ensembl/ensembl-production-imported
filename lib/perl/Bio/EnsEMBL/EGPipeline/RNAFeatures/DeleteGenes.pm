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

Bio::EnsEMBL::EGPipeline::RNAFeatures::DeleteGenes

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::RNAFeatures::DeleteGenes;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub run {
  my ($self) = @_;
  my $logic_name = $self->param_required('logic_name');
  
  my $dba  = $self->core_dba();
  my $ga   = $dba->get_adaptor('Gene');
  my $dbea = $dba->get_adaptor('DBEntry');
  my $aa   = $dba->get_adaptor('Attribute');
  
  my @genes = @{ $ga->fetch_all_by_logic_name($logic_name) };
  
  # Unfortunately can't just do $ga->remove, because that also
  # deletes all of the dna_align_features...
  foreach my $gene (@genes) {
    $self->remove_gene($dba, $dbea, $aa, $gene);
  }
}

sub remove_gene {
  my ($self, $dba, $dbea, $aa, $gene) = @_;
  
  # Remove object xrefs
  foreach my $dbe (@{ $gene->get_all_DBEntries() }) {
    $dbea->remove_from_object($dbe, $gene, 'Gene');
  }
  
  # Remove attributes
  $aa->remove_from_Gene($gene);
  
  # Remove transcripts
  foreach my $transcript (@{ $gene->get_all_Transcripts() }) {
    $self->remove_transcript($dba, $dbea, $aa, $transcript);
  }
  
  my $sth = $dba->dbc->prepare("DELETE FROM gene WHERE gene_id = ? ");
  $sth->execute($gene->dbID);
}

sub remove_transcript {
  my ($self, $dba, $dbea, $aa, $transcript) = @_;
  
  # Remove object xrefs
  foreach my $dbe (@{ $transcript->get_all_DBEntries() }) {
    $dbea->remove_from_object($dbe, $transcript, 'Transcript');
  }
  
  # Remove attributes
  $aa->remove_from_Transcript($transcript);
  
  # Remove the translation if it exists
  if (defined $transcript->translation) {
    my $ta = $dba->get_adaptor('Translation');
    $ta->remove($transcript->translation);
  }
  
  # Remove exons
  foreach my $exon (@{ $transcript->get_all_Exons() }) {
    # Only remove the exon if this is the last transcript to reference it
    
    my $sth = $dba->dbc->prepare("SELECT count(*) FROM exon_transcript WHERE exon_id = ?" );
    $sth->execute($exon->dbID);
    my ($count) = $sth->fetchrow_array();
    
    if ($count == 1) {
      $self->remove_exon($dba, $exon);
    }
  }
  
  # Remove exon/transcript links
  my $sth = $dba->dbc->prepare("DELETE FROM exon_transcript WHERE transcript_id = ?");
  $sth->execute($transcript->dbID);

  # Remove supporting feature links
  $sth = $dba->dbc->prepare("DELETE FROM transcript_supporting_feature WHERE transcript_id = ?");
  $sth->execute($transcript->dbID);

  # Remove transcript
  $sth = $dba->dbc->prepare( "DELETE FROM transcript WHERE transcript_id = ?" );
  $sth->execute($transcript->dbID);
}

sub remove_exon {
  my ($self, $dba, $exon) = @_;
  
  # Remove supporting feature links
  my $sth = $dba->dbc->prepare("DELETE FROM supporting_feature WHERE exon_id = ?");
  $sth->execute($exon->dbID);
  
  # Remove exon
  $sth = $dba->dbc->prepare("DELETE FROM exon WHERE exon_id = ?");
  $sth->execute($exon->dbID);
}

1;
