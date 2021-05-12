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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::ConvertSpace;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Capture::Tiny ':all';
use Path::Tiny qw(path);
use File::Basename qw(dirname);
use File::Spec::Functions qw(catdir);

sub param_defaults {
  my ($self) = @_;
  
  return {
    'threads' => 1,
  };
}

sub run {
  my ($self) = @_;
  
  my $seq1 = $self->param_required('seq_file_1');
  my $seq2 = $self->param('seq_file_2');

  my $space = space_check($seq1);
  
  if ($space eq 'color') {
    # Convert
    convert_to_base_space($seq1);
    convert_to_base_space($seq2) if $seq2;
  }
  
  my $output = {
    seq_file_1 => $seq1,
    seq_file_2 => $seq2,
  };
  $self->dataflow_output_id($output,  2);
}

sub space_check {
  my ($path) = @_;
  
  open my $fh, '<:gzip', $path;
  my $id_line = readline $fh;
  my $sequence_line = readline $fh;
  close $fh;
  
  return sequence_space($sequence_line);
}

sub sequence_space {
  my ($seq) = @_;

  if ($seq =~ /^[CGTAN]+$/) {
    return "base";
  } elsif ($seq =~ /^[CGTA][0-3\.]+$/) {
    return "color";
  } else {
    die "Could not determine the color space: '$seq'";
  }
}

sub convert_to_base_space {
  my ($inpath) = @_;
  return if not $inpath;
  
  my $outpath = $inpath . ".base";
  
  warn("Convert $inpath to $outpath");
  open my $infh, '<:gzip', $inpath;
  open my $outfh, '>:gzip', $outpath;
  
  while (my $id_line1 = readline $infh) {
    # Get the 4 lines of the read
    my $color_sequence = readline $infh;
    my $id_line2 = readline $infh;
    my $color_quality = readline $infh;
    
    my $base_sequence = convert_sequence($color_sequence);
    my $base_quality = convert_quality($color_quality);
    
    print $outfh $id_line1;
    print $outfh $base_sequence;
    print $outfh $id_line2;
    print $outfh $base_quality;
  }
  close $infh;
  
  rename $outpath, $inpath;
}

sub convert_sequence {
  my ($cseq) = @_;
  
  chomp $cseq;
  my $bseq;
  
  my %base_map = (
    A => [qw(A C G T)],
    C => [qw(C A T G)],
    G => [qw(G T A C)],
    T => [qw(T G C A)],
  );
  
  my ($first, @colors) = split //, $cseq;
  my $map = $base_map{$first};
  
  my $fail = 0;
  my $nucl = $first;
  for my $c (@colors) {
    if ($c eq '.' or int($c) > 3) {
      $fail = 1;
    }
    
    if ($fail) {
      $bseq .= "N";
    } else {
      $nucl = $base_map{$nucl}[$c];
      $bseq .= $nucl;
    }
  }
  
  return $bseq . "\n";
}

sub convert_quality {
  my ($cqual) = @_;
  
  my $bqual = $cqual;
  $bqual =~ s/^.//;
  
  return $bqual;
}

1;
