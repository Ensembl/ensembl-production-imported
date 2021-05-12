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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::SplitUnique;

use strict;
use warnings;
use Capture::Tiny ':all';

use base ('Bio::EnsEMBL::EGPipeline::BRC4Aligner::Base');

sub run {
  my ($self) = @_;

  $self->dbc and $self->dbc->disconnect_if_idle();
  
  my $bam_file = $self->param_required('bam_file');

  my @bam_files;

  # Split in two files: unique maps, and multiple maps
  my @cmds;
  my @bams;
  for my $number (qw{ unique non_unique }) {
    my ($bam_output, $cmd) = $self->split_unique($bam_file, $number);
    index_bam($bam_output);
    push @bams, $bam_output;
    push @cmds, $cmd;
  }
  $self->param("cmds", \@cmds);
  $self->param("bams", \@bams);
}

sub write_output {
  my ($self) = @_;
  
  my $cmds = $self->param("cmds");
  my $bams = $self->param("bams");

  for my $cmd (@$cmds) {
    my $align_cmds = {
      cmds => $cmd,
      sample_name => $self->param('sample_name'),
    };
    $self->store_align_cmds($align_cmds);
  }

  for my $bam (@$bams) {
    $self->dataflow_output_id({ bam_file => $bam },  2);
  }
}

sub split_unique {
  my ($self, $bam_input, $number) = @_;

  my $bam_output = $bam_input;
  if ($bam_input =~ /results\.bam/) {
    $bam_output =~ s/results\.bam/${number}_results.bam/;
  } else {
    $bam_output =~ s/\.bam/_${number}.bam/;
  }

  # Prepare the command for the RSeQ script
  my $grep_command;

  if ($number eq 'unique') {
    #$grep_command = "grep -P '^\@|NH:i:1(\\s|\$)'";
    $grep_command = "grep -P '^\@|NH:i:1\\b'";
  } elsif ($number eq 'non_unique') {
    $grep_command = "grep -Pv 'NH:i:1\\b'";
  } else {
    die "Unrecognized uniqueness to split: $number";
  }

  my $cmd = "samtools view -h $bam_input | $grep_command | samtools view -h -b > $bam_output";

  my ($stdout, $stderr, $exit) = capture {
    system( $cmd );
  };

  if ($exit != 0) {
    die("unique split failed ($exit): '$stderr'");
  }

  return ($bam_output, $cmd);
}

sub index_bam {
  my ($bam) = @_;

  my $cmd = "samtools index $bam";
  my ($stdout, $stderr, $exit) = capture {
    system($cmd);
  };

  if ($exit != 0) {
    die("Index failed: $stderr");
  }

  return $cmd;
}

1;
