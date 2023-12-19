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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::CalculateTPM;
use strict;
use warnings;
use autodie;
use base ('Bio::EnsEMBL::EGPipeline::BRC4Aligner::Base');


sub run {
  my ($self) = @_;

  my $sum_htseq_file = $self->param_required('sum_file');
  my $fl_file = $self->param_required('fl_file');
  my $tpm_htseq_file = $sum_htseq_file . ".tpm";
  my $rpk_htseq_file = $sum_htseq_file . ".rpk";

  my @rpk_args = ( "bash", "-c", "join <(sort -k1,1 $sum_htseq_file) <(sort -k1,1 $fl_file) | awk '{print \$1\"\t\"\$2/(\$3/1000)}' > $rpk_htseq_file" );
  system(@rpk_args);

  open my $fh, '<', $rpk_htseq_file;

  my $total = 0;

  while (<$fh>) {
    my @data = split /\t/, $_;

    die "The $rpk_htseq_file contains more than 3 columns. Exiting." unless scalar(@data) == 2;

    my $gene_id     	  = $data[0];
    my $gene_factor   	  = $data[1];

    $total += $gene_factor;
  }

  my $million_factor = $total/1000000;

  my @tpm_args = ( "bash", "-c", "awk -v million_factor=\"$million_factor\" '{print \$1\"\t\"\$2/million_factor}' $rpk_htseq_file > $tpm_htseq_file");
  system(@tpm_args);

  unlink $rpk_htseq_file;

  $self->param("tpm_htseq_file", $tpm_htseq_file);
}

sub write_output {

  my ($self) = @_;

  my $tpm_htseq_file = $self->param("tpm_htseq_file");

  $self->dataflow_output_id({ tpm_htseq_file => $tpm_htseq_file }, 2);

}

1;
