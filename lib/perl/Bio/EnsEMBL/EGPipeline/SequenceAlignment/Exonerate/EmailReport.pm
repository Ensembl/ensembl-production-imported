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

Bio::EnsEMBL::EGPipeline::SequenceAlignment::Exonerate::EmailReport

=head1 DESCRIPTION

Run a few useful queries on the results of exonerate, for a given species.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::SequenceAlignment::Exonerate::EmailReport;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EmailReport');

use Bio::SeqIO;

sub param_defaults {
  my $self = shift @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'db_type'  => 'otherfeatures',
  };
}

sub fetch_input {
  my ($self) = @_;
  
  my $species          = $self->param_required('species');
  my $logic_name       = $self->param_required('logic_name');
  my $seq_file         = $self->param('seq_file');
  my $seq_file_species = $self->param('seq_file_species');
  my $data_type        = $self->param('data_type');
  my $coverage         = $self->param('coverage');
  my $percent_id       = $self->param('percent_id');
  
  if (exists $$seq_file_species{$species}) {
    push @$seq_file, $$seq_file_species{$species};
  }
  my $fasta_files = join(', ', @$seq_file);
  
  my $seq_type = $data_type eq 'protein' ? 'protein' : 'dna';
  
  my ($seq_count, $seq_length) = (0, 0);
  foreach my $fasta_file (@$seq_file) {
    my $seqs = Bio::SeqIO->new(-format => 'Fasta', -file => $fasta_file);
    while (my $seq = $seqs->next_seq) {
      $seq_count++;
      $seq_length += $seq->length;
    }
  }
  my $length = sprintf("%.0f", $seq_length/1000000);
  my $mean_length = sprintf("%.1f", $seq_length/($seq_count*1000));
  my $length_unit = $data_type eq 'protein' ? 'million amino acids' : 'Mb';
  my $mean_length_unit = $data_type eq 'protein' ? 'amino acids' : 'Kb';
  
  my $dba = $self->get_DBAdaptor($self->param('db_type'));
  my $dbh = $dba->dbc->db_handle();
  
  my $sql =
    'SELECT COUNT(distinct hit_name) FROM '.
    $seq_type.'_align_feature INNER JOIN analysis USING (analysis_id) '.
    'WHERE logic_name = "'.$logic_name.'";';
  my ($unique_hits) = $dbh->selectrow_array($sql);
  
  my $seq_hit_pcage = sprintf("%.0f", ($unique_hits/$seq_count)*100);
  
  my $text = 
    "The exonerate pipeline has completed for $species, ".
    "using the sequence file(s): $fasta_files. ".
    "That file has $seq_count sequences and a total length of $length $length_unit, ".
    "giving an average sequence length of $mean_length $mean_length_unit. ".
    "Of that total, $unique_hits ($seq_hit_pcage%) were mapped to the genome.\n\n";
  
  if ($self->param('make_genes')) {
    $sql =
    'SELECT '.
    'COUNT(distinct gene_id) AS genes, '.
    'COUNT(distinct transcript_id) AS transcripts, '.
    'COUNT(exon_id) AS exons FROM '.
    'transcript INNER JOIN '.
    'exon_transcript USING (transcript_id) INNER JOIN '.
    'analysis USING (analysis_id) '.
    'WHERE logic_name = "'.$logic_name.'";';
    my ($genes, $transcripts, $exons) = $dbh->selectrow_array($sql);
    
    $text .=
      "With thresholds of $coverage% coverage and $percent_id% sequence ".
      "identity, the pipeline generated $genes genes, ".
      "$transcripts transcripts, and $exons exons.\n\n";
  }
  
  $text .=
    "Fond regards,\nThe Exonerate Pipeline\n";
  
  $self->param('text', $text);
}

1;
