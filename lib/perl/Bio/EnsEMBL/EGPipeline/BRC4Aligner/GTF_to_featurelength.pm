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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::GTF_to_featurelength;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub run {
  my ($self) = @_;

  my $gtf_file = $self->param_required('gtf_file');
  my $fl_file = $self->param_required('fl_file');

  my %el_dict;

  open my $gtf, "<", $gtf_file or die("$gtf_file: $!");
  open my $el, ">", $fl_file or die("$fl_file: $!");

  while (my $line = readline($gtf)) {
      next if $line =~ /^#/;
      chomp $line;

      my ($chr, $source, $type, $start, $end, $score, $strand, $phase, $attr_list) = split /\t/, $line;

      # We only need exons positions here
      next if $type ne 'exon';

      # Get exon's length
      my $exon_length = abs($start - $end);

      # The exon id is needed
      my %attr = get_attributes($attr_list);
      my $gene_id = $attr{"gene_id"};
      $gene_id =~ s/"//g;

      # Increase the genes length value by the exon's length in the dict:
      $el_dict{$gene_id} += $exon_length;
  }

  my $gene_id;
  my $exons_total_length;
while ( ($gene_id,$exons_total_length) = each %el_dict ) {

    my @out_line = ($gene_id, $exons_total_length);

    print $el join("\t", @out_line) . "\n";
  }
  close $gtf;
  close $el;

  $self->param('fl_file', $fl_file);
}

sub write_output {
  my ($self) = @_;

  $self->dataflow_output_id({ fl_file => $self->param('fl_file') },  2);
}

sub get_attributes {
  my ($attr_list) = @_;

  my @attrs = split /; /, $attr_list;

  my %attrs_dict;
  for my $attr (@attrs) {
    if ($attr =~ /^(.+) (.+)$/) {
      $attrs_dict{$1} = $2;
    }
  }
  return %attrs_dict;
}

1;
