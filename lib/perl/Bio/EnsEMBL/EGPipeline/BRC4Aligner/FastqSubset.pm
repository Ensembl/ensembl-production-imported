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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::FastqSubset;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');


sub param_defaults {
  my ($self) = @_;
  
  return {
    n_samples => 1_000_000
  };
}

sub run {
  my ($self) = @_;
  
  my $seq1 = $self->param_required('seq_file_1');
  my $seq2 = $self->param('seq_file_2');
  my $n_samples = $self->param_required('n_samples');

  my $sub_seq1 = $seq1; $sub_seq1 =~ s/(\..+$)/_subset$1/;
  $self->subset_sequences($seq1, $sub_seq1, $n_samples);
  $self->param('subset_seq_file_1', $sub_seq1);

  if ($seq2) {
    my $sub_seq2 = $seq2; $sub_seq2 =~ s/(\..+$)/_subset$1/;
    $self->subset_sequences($seq2, $sub_seq2, $n_samples);
    $self->param('subset_seq_file_2', $sub_seq2);
  } else {
    $self->param('subset_seq_file_2', $seq2);
  }
}

sub write_output {
  my ($self) = @_;
  
  my $output = { subset_seq_file_1 => $self->param('subset_seq_file_1') };

  if (defined $self->param('subset_seq_file_2')) {
    $output->{subset_seq_file_2} = $self->param('subset_seq_file_2');
  } else {
    $output->{subset_seq_file_2} = undef;
  }
  my $sam_file = $self->param('subset_seq_file_1') . ".sam";
  $output->{sam_file} = $sam_file;
  $self->dataflow_output_id($output,  1);
}

sub subset_sequences {
  my ($self, $seq_file, $sample_file, $n_samples) = @_;

  open my $in, "<:gzip", $seq_file or die ("Can't open $seq_file to read: $!");
  open my $out, ">:gzip", $sample_file or die ("Can't open $sample_file to write: $!");
  
  my $read_count = 0;
  my $line_count = 0;

  while (my $line = readline $in) {

    # One read
    if ($line_count++ % 4 == 0) {
      $read_count++;
      last if $read_count > $n_samples; 
    }
    print $out $line;
  }
  close $in;
  close $out;
}

1;
