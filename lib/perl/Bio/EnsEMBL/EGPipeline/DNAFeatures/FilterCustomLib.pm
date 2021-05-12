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

package Bio::EnsEMBL::EGPipeline::DNAFeatures::FilterCustomLib;

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
    'threads' => 1,
  };
  
}

sub run {
  my ($self) = @_;
  
  my $transcripts_file = $self->param_required('transcripts_file');
  my $repeat_lib_to_filter = $self->param_required('lib_to_filter');
  my $filtered_repeat_lib = $self->param_required('filtered_lib');
  
  my $rep_hits = $self->hit_repeat_lib($repeat_lib_to_filter, $transcripts_file);
  $self->remove_hits($repeat_lib_to_filter, $filtered_repeat_lib, $rep_hits);
}

sub hit_repeat_lib {
  my ($self, $lib_to_filter, $transcripts_file) = @_;
  
  my $threads = $self->param('threads');
  
  # Create a blast DB from the repeat library
  my $blast_db = $transcripts_file . ".db";
  my $db_cmd = "makeblastdb -in $transcripts_file -dbtype nucl -input_type fasta -out $blast_db";
  print STDERR "Blast db command: $db_cmd\n";
  system($db_cmd);
  
  # Blast the proteins against the repeat lib
  my $blast_result = $blast_db . ".lib_vs_cleantr";
  my $blast_cmd = "blastn -db $blast_db -query $lib_to_filter";
  my @parameters = (
    '-evalue 1e-10',
    '-culling_limit 2',
    '-max_target_seqs 25',
    '-outfmt "6 qseqid bitscore std"',
    "-num_threads $threads",
    "-out $blast_result"
  );
  $blast_cmd .= " " . join(" ", @parameters);
  print STDERR "Blast command: $blast_cmd\n";
  system($blast_cmd);

  my $hits = $self->get_blast_hits($blast_result);
  print STDERR scalar(keys %$hits) . " hits\n";
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
  
  return \%hits;
}

sub remove_hits {
  my ($self, $input, $output, $hits) = @_;
  
  my $new_file = $output . ".filter";
  my $seq_in = Bio::SeqIO->new(-file => $input, -format => "fasta");
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
  
  print STDERR "$skipped sequences were removed\n";
  
  rename $new_file, $output;
}

1;
