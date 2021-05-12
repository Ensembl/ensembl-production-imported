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

package Bio::EnsEMBL::EGPipeline::DNAFeatures::RepbaseTranslationsToIgnore;

use strict;
use warnings;

use File::Spec::Functions qw(catfile catdir);
use File::Basename qw(basename dirname);
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::SeqIO;

sub param_defaults {
  my ($self) = @_;
  
  return {
    %{$self->SUPER::param_defaults},
    ignore_biotypes => ['transposable_element'],
    'threads' => 1,
  };
  
}

sub run {
  my ($self) = @_;
  
  my $species_repeat_lib = $self->param_required('species_repeat_lib');

  my $transcripts_file = $self->param_required('transcripts_file');
  my $translations_file = $self->param_required('translations_file');
  
  my $translation_hits = $self->hit_translations($translations_file, $species_repeat_lib);
  $self->filter_transcripts($transcripts_file, $translation_hits);
}

sub hit_translations {
  my ($self, $fasta_input, $rep_lib) = @_;
  
  my $threads = $self->param('threads');
  
  # Create a blast DB from the repeat library
  my $blast_db = $rep_lib . ".db";
  my $db_cmd = "makeblastdb -in $rep_lib -dbtype nucl -input_type fasta -out $blast_db";
  print STDERR "Blast db command: $db_cmd\n";
  system($db_cmd);
  
  # Blast the proteins against the repeat lib
  my $blast_result = $blast_db . ".pep_vs_repbase";
  my $blast_cmd = "tblastn -db $blast_db -query $fasta_input";
  my @parameters = (
    '-evalue 1e-5',
    '-culling_limit 2',
    '-max_target_seqs 10',
    '-outfmt "6 qseqid bitscore std"',
    "-num_threads $threads",
    "-out $blast_result"
  );
  $blast_cmd .= " " . join(" ", @parameters);
  print STDERR "Blast command: $blast_cmd\n";
  system($blast_cmd);

  my $hits = $self->get_blast_hits($blast_result);
  print STDERR scalar(@$hits) . " hits\n";
  return $hits;
}

sub get_blast_hits {
  my ($self, $blast_result) = @_;
  
  my %hits;
  open my $blast_fh, "<", $blast_result;
  while (my $line = readline $blast_fh) {
    my ($id) = split /\s+/, $line;
    $hits{$id} = 1;
  }
  close $blast_fh;
  
  return [keys %hits];
}

sub filter_transcripts {
  my ($self, $transcripts_file, $hits) = @_;
  
  # Map the translation ids to the transcripts
  my $tr_map = $self->get_translations_map();
  my %mapped_hits = map { $_ => 1 } map { exists $tr_map->{$_} ? $tr_map->{$_} : die "No transcript for translation $_" } @$hits;

  # Also remove any ignored biotypes
  my $to_ignore = $self->get_transcripts_to_ignore();
  
  my %to_remove = (%mapped_hits, %$to_ignore);
  
  # Remove hit transcripts
  $self->remove_hit_transcripts($transcripts_file, \%to_remove);

}

sub get_translations_map {
  my ($self) = @_;
  
  print STDERR "Get translation -> transcript map\n";
  
  my $dba = $self->get_DBAdaptor('core');
  my $tla = $dba->get_adaptor('translation');

  my %map;
  my $translations = $tla->fetch_all();
  print STDERR scalar(@$translations) . " translations\n";
  for my $translation (@$translations) {
    my $transcript = $translation->transcript;
    $map{$translation->stable_id} = $transcript->stable_id;
  }
  $dba->dbc->disconnect_if_idle();
  
  return \%map;
}

sub get_transcripts_to_ignore {
  my ($self) = @_;
  
  my $biotypes_to_ignore = $self->param("ignore_biotypes");
  return if not $biotypes_to_ignore;
  
  print STDERR "Get list of transcripts to ignore (based on biotype)\n";
  
  my $dba = $self->get_DBAdaptor('core');
  my $tra = $dba->get_adaptor('transcript');

  my $transcripts = $tra->fetch_all_by_biotype($biotypes_to_ignore);
  print STDERR scalar(@$transcripts) . " transcripts to ignore\n";
  my %transcript_ids = map { $_->stable_id => 1 } @$transcripts;
  
  return \%transcript_ids;
}

sub remove_hit_transcripts {
  my ($self, $transcripts_file, $hits) = @_;
  
  my $new_file = $transcripts_file . ".filter";
  my $old_file = $transcripts_file . ".prefilter";
  my $seq_in = Bio::SeqIO->new(-file => $transcripts_file, -format => "fasta");
  my $seq_out = Bio::SeqIO->new(-file => ">$new_file", -format => "fasta");
  
  my $skipped = 0;
  while (my $inseq = $seq_in->next_seq) {
    if (exists $hits->{$inseq->id}) {
      $skipped++;
      next;
    } else {
      $seq_out->write_seq($inseq);
    }
  }
  
  print STDERR "$skipped transcripts were removed\n";
  
  rename $transcripts_file, $old_file;
  rename $new_file, $transcripts_file;
}

1;
