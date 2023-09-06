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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::SplitReadStrand;

use strict;
use warnings;
use Capture::Tiny ':all';

use base ('Bio::EnsEMBL::EGPipeline::BRC4Aligner::Base');

sub run {
  my ($self) = @_;

  $self->dbc and $self->dbc->disconnect_if_idle();
  
  my $bam_file = $self->param_required('bam_file');
  my $aligner_metadata = $self->param_required('aligner_metadata');
  my $is_paired = $aligner_metadata->{'is_paired'} // $self->param{'is_paired'};
  my $is_stranded = $aligner_metadata->{'is_stranded'};
  my $strand_direction = $aligner_metadata->{'strand_direction'};

  # Split in two strand files
  if ($is_stranded) {
    die("Data is stranded but has no direction?") if not $strand_direction;
    
    my $cmds;
    if ($is_paired) {
      if ($strand_direction eq 'forward') {
        $cmds = $self->split_paired_forward($bam_file);
      } elsif ($strand_direction eq 'reverse') {
        $cmds = $self->split_paired_reverse($bam_file);
      } else {
        die("Strand direction needed: $strand_direction");
      }
    } else {
      if ($strand_direction eq 'forward') {
        $cmds = $self->split_single_forward($bam_file);
      } elsif ($strand_direction eq 'reverse') {
        $cmds = $self->split_single_reverse($bam_file);
      }
    }

    # Store align commands
      my $align_cmds = {
        cmds => $cmds,
        sample_name => $self->param('sample_name'),
      };
      $self->store_align_cmds($align_cmds);

  } else {
    print("The reads are not strand-specific: skip splitting");
    $self->dataflow_output_id({ bam_file => $bam_file },  2);
  }
}

sub split_paired_forward {
  my ($self, $bam_input) = @_;

  # Prepare files
  my $fwd = $bam_input;
  $fwd =~ s/\.bam/.firststrand.bam/;
  my $fwd1 = $bam_input;
  my $fwd2 = $bam_input;
  $fwd1 =~ s/\.bam/_fwd1.bam/;
  $fwd2 =~ s/\.bam/_fwd2.bam/;

  my $rev = $bam_input;
  $rev =~ s/\.bam/.secondstrand.bam/;
  my $rev1 = $bam_input;
  my $rev2 = $bam_input;
  $rev1 =~ s/\.bam/_rev1.bam/;
  $rev2 =~ s/\.bam/_rev2.bam/;

  my @cmds;

  # Forward
  push @cmds, extract_bam($bam_input, $fwd1, "-f 64  -F 20");
  push @cmds, extract_bam($bam_input, $fwd2, "-f 144 -F 4");
  push @cmds, merge_bams($fwd1, $fwd2, $fwd);
  push @cmds, index_bam($fwd);
  $self->dataflow_output_id({ bam_file => $fwd, strand => 'firststrand' },  2);

  # Reverse
  push @cmds, extract_bam($bam_input, $rev1, "-f 80  -F 4");
  push @cmds, extract_bam($bam_input, $rev2, "-f 128 -F 20");
  push @cmds, merge_bams($rev1, $rev2, $rev);
  push @cmds, index_bam($rev);
  $self->dataflow_output_id({ bam_file => $rev, strand => 'secondstrand' },  2);

  # Clean up
  my @to_delete = ($fwd1, $fwd2, $rev1, $rev2);
  for my $to_del (@to_delete) {
    unlink $to_del;
  }

  return \@cmds;
}

sub split_paired_reverse {
  my ($self, $bam_input) = @_;

  # Prepare files
  my $fwd = $bam_input;
  $fwd =~ s/\.bam/.firststrand.bam/;
  my $fwd1 = $bam_input;
  my $fwd2 = $bam_input;
  $fwd1 =~ s/\.bam/_fwd1.bam/;
  $fwd2 =~ s/\.bam/_fwd2.bam/;

  my $rev = $bam_input;
  $rev =~ s/\.bam/.secondstrand.bam/;
  my $rev1 = $bam_input;
  my $rev2 = $bam_input;
  $rev1 =~ s/\.bam/_rev1.bam/;
  $rev2 =~ s/\.bam/_rev2.bam/;

  my @cmds;

  # Forward
  push @cmds, extract_bam($bam_input, $fwd1, "-f 80  -F 4");
  push @cmds, extract_bam($bam_input, $fwd2, "-f 128 -F 20");
  push @cmds, merge_bams($fwd1, $fwd2, $fwd);
  push @cmds, index_bam($fwd);
  $self->dataflow_output_id({ bam_file => $fwd, strand => 'firststrand' },  2);

  # Reverse
  push @cmds, extract_bam($bam_input, $rev1, "-f 64  -F 20");
  push @cmds, extract_bam($bam_input, $rev2, "-f 144 -F 4");
  push @cmds, merge_bams($rev1, $rev2, $rev);
  push @cmds, index_bam($rev);
  $self->dataflow_output_id({ bam_file => $rev, strand => 'secondstrand' },  2);

  # Clean up
  my @to_delete = ($fwd1, $fwd2, $rev1, $rev2);
  for my $to_del (@to_delete) {
    unlink $to_del;
  }

  return \@cmds;
}

sub split_single_forward {
  my ($self, $bam_input) = @_;

  # Prepare files
  my $fwd = $bam_input;
  $fwd =~ s/\.bam/.firststrand.bam/;

  my $rev = $bam_input;
  $rev =~ s/\.bam/.secondstrand.bam/;

  my @cmds;

  # Forward
  push @cmds, extract_bam($bam_input, $fwd, "-F 20");
  push @cmds, index_bam($fwd);
  $self->dataflow_output_id({ bam_file => $fwd, strand => 'firststrand' },  2);

  # Reverse
  push @cmds, extract_bam($bam_input, $rev, "-f 16 -F 4");
  push @cmds, index_bam($rev);
  $self->dataflow_output_id({ bam_file => $rev, strand => 'secondstrand' },  2);

  return \@cmds;
}

sub split_single_reverse {
  my ($self, $bam_input) = @_;

  # Prepare files
  my $fwd = $bam_input;
  $fwd =~ s/\.bam/.firststrand.bam/;

  my $rev = $bam_input;
  $rev =~ s/\.bam/.secondstrand.bam/;

  my @cmds;

  # Forward
  push @cmds, extract_bam($bam_input, $fwd, "-f 16 -F 4");
  push @cmds, index_bam($fwd);
  $self->dataflow_output_id({ bam_file => $fwd, strand => 'firststrand' },  2);

  # Reverse
  push @cmds, extract_bam($bam_input, $rev, "-F 20");
  push @cmds, index_bam($rev);
  $self->dataflow_output_id({ bam_file => $rev, strand => 'secondstrand' },  2);

  return \@cmds;
}

sub extract_bam {
  my ($bam_input, $bam_output, $filters) = @_;

  my $cmd = "samtools view -h -b $filters $bam_input > $bam_output";

  my ($stdout, $stderr, $exit) = capture {
    system( $cmd );
  };

  if ($exit != 0) {
    die("Forward split failed: $stderr");
  }

  return $cmd;
}

sub merge_bams {
  my ($bam1, $bam2, $bam_out) = @_;

  my $cmd = "samtools merge -f $bam_out $bam1 $bam2";

  my ($stdout, $stderr, $exit) = capture {
    system( $cmd );
  };

  if ($exit != 0) {
    die("Forward split failed: $stderr");
  }

  return $cmd;
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
