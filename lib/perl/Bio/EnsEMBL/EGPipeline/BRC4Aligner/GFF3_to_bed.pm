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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::GFF3_to_bed;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub run {
  my ($self) = @_;
  
  my $gff_file = $self->param_required('gff_file');
  my $bed_file = $self->param_required('bed_file');
  
  open my $gff, "<", $gff_file or die("$gff_file: $!");
  open my $bed, ">", $bed_file or die("$bed_file: $!");
  print $bed 'track name=exons description="Exons"'."\n";

  while (my $line = readline($gff)) {
    next if $line =~ /^#/;
    chomp $line;

    my ($chr, $source, $type, $start, $end, $score, $strand, $phase, $attr_list) = split /\t/, $line;

    # We only need exons positions here
    next if $type ne 'exon';

    # The exon id is needed
    my %attr = get_attributes($attr_list);
    my $id = $attr{ID} || $attr{Name};
    $score = 0;
    my @bed_line = ($chr, $start, $end, $id, $score, $strand);

    print $bed join("\t", @bed_line) . "\n";
  }
  close $gff;
  close $bed;
  
  $self->param('bed_file', $bed_file);
}

sub write_output {
  my ($self) = @_;
  
  $self->dataflow_output_id({ bed_file => $self->param('bed_file') },  2);
}

sub get_attributes {
  my ($attr_list) = @_;

  my @attrs = split /;/, $attr_list;

  my %attrs_dict;
  for my $attr (@attrs) {
    if ($attr =~ /^(.+)=(.+)$/) {
      $attrs_dict{$1} = $2;
    }
  }
  return %attrs_dict;
}

1;
