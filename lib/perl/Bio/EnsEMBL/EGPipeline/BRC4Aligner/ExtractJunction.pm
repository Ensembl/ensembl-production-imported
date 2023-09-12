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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::ExtractJunction;

use strict;
use warnings;
use v5.14;
use autodie;

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');
use File::Spec::Functions qw(catdir);

sub run {
  my ($self) = @_;
  
  my $bam_file = $self->param_required('bam_file');
  my $aligner_metadata = $self->param_required('aligner_metadata');

  # Junction file
  my $bam_results_dir = $self->param_required('results_dir');
  my $junction_file = catdir($bam_results_dir, 'junctions.tab');

  # Get the splice junctions
  my $min_size = 0;
  $self->get_splice_junctions($bam_file, $junction_file, $min_size);
}

sub get_splice_junctions {
  my ($self, $bam, $junction_file, $min_size) = @_;

  my %junctions;
  open INBAM, "samtools view -F4 $bam|";
  while(my $line = readline INBAM) {
    chomp $line;

    my ($name, $flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual, @tags_list) = split /\t/, $line;

    # We only want reads with splice junctions
    next unless $cigar =~ /N/;

    # Check tags
    my %tag = $self->get_tags(\@tags_list);
    next unless $tag{XS};
    die "Tag not found: NH in $line" unless $tag{NH};
    
    # Mapper type?
    my $mapper = "nu";
    if ($tag{NH} == 1) {
      $mapper = "unique";
    }

    # Get splice-junction coordinates and size
    my $running_offset = 0;
    while($cigar =~ s/^(\d+)([^\d])//) {
      my $num = $1;
      my $type = $2;

      if($type eq 'N') {
        my $start = $pos + $running_offset;
        my $end = $start + $num - 1;
        my $junction = "$rname:$start-$end";

        if($num >= $min_size) {
          $junctions{$junction}->{ $tag{XS} }->{$mapper}++;
        }
      }

      if($type eq 'N' || $type eq 'M' || $type eq 'D') {
        $running_offset = $running_offset + $num;
      }
    }
  }
  close INBAM;

  # Print the junctions
  
  open OUTJUNC, ">", $junction_file;
  my @fields = qw(Junction Strand Unique  NU);
  print OUTJUNC join("\t", @fields) . "\n";

  for my $junction (sort keys %junctions) {
    foreach my $strand (sort keys %{$junctions{$junction}}) {
      my $junc = $junctions{$junction}->{$strand};

      my $unique_count = $junc->{unique} // 0;
      my $nu_count = $junc->{nu} // 0;

      my @data = ($junction, $strand, $unique_count, $nu_count);
      print OUTJUNC join("\t", @data) . "\n";

    }
  }
}

sub get_tags {
  my ($self, $tags_list) = @_;

  my %tags;
  for my $tag (@$tags_list) {
    if ($tag =~ /^(..):.:(.+)$/) {
      $tags{$1} = $2;
    }
  }

  return %tags;
}

1;
