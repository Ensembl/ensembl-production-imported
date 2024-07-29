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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::BamStats;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Capture::Tiny ':all';
use Path::Tiny qw(path);
use File::Basename qw(dirname basename);
use File::Spec::Functions qw(catdir);

sub write_output {
  my ($self) = @_;

  my $bam = $self->param_required('bam_file');
  my $res_dir = $self->param('results_dir');
  my $length_file = $self->param_required('length_file');
  my $global_total_reads = $self->param('total_reads');

  # Run
  my $stats = $self->run_bam_stats($bam, $length_file, $global_total_reads);

  use Data::Dumper;
  print(Dumper($stats)."\n");
  
  # Only continue if reads have been aligned
  if ($stats->{number_reads_mapped} > 0) {
    my $data = { bam_stats => $stats };

    # No total reads yet, means that this is the main bam file, need to propagate total reads
    if (not $global_total_reads) {
      $data->{total_reads} = $stats->{total_reads};
    }

    $self->dataflow_output_id($data, 2);
  } else {
    if ($res_dir) {
      my $log_file = catdir($res_dir, 'log.txt');
      print "WARNING: No reads aligned";
      open my $LOG, ">>", $log_file;
      close $LOG;
    }
  }

  # Store bam_stats in all cases
  $self->dataflow_output_id({ bam_stats => $stats }, 3);
}

sub run_bam_stats {
  my ($self, $bam, $length_file, $global_total_reads) = @_;

  # Get coverage
  my $coverage = $self->get_coverage($bam, $length_file);

  # Get stats from sam
  my $stats = $self->get_samtools_stats($bam);
  my $total_reads = $stats->{'raw total sequences'};
  my $number_mapped = $stats->{'reads mapped'};
  my $average_read_length = $stats->{'average length'};
  my $percent_mapped = $global_total_reads ? $number_mapped / $global_total_reads : ($total_reads ? $number_mapped / $total_reads : 0);

  my %stats = (
    file => basename($bam),
    total_reads => $total_reads,
    coverage => $coverage,
    mapped => $percent_mapped,
    number_reads_mapped => $number_mapped,
    average_read_length => $average_read_length,
  );
  my $aligner_metadata = $self->param_required('aligner_metadata');
  $stats{number_pairs_mapped} = $stats->{'reads properly paired'} if $aligner_metadata->{'is_paired'};

  return \%stats;
}

sub get_coverage {
  my ($self, $bam, $length_file) = @_;

  my $cmd = "bedtools genomecov -ibam $bam -g $length_file";
  my ($stdout, $stderr, $exit) = capture {
    system($cmd);
  };

  if ($exit != 0) {
    die("Index failed: $stderr, $stdout");
  }

  my $genomeCoverage = 0;
  my $count = 0;
  for my $line (split /\n/, $stdout) {
    if ($line =~ /^genome/) {
      my ($identifier, $depth, $freq, $size, $proportion) = split(/\t/, $line);
      $genomeCoverage += ($depth * $freq);
      $count += $freq;
    }
  }

  return $count ? $genomeCoverage / $count : 0;
}

sub get_samtools_stats {
  my ($self, $bam) = @_;

  my $cmd = "samtools stats $bam | grep ^SN | cut -f 2-";
  my ($stdout, $stderr, $exit) = capture {
    system($cmd);
  };

  if ($exit != 0) {
    die("Index failed: $stderr");
  }

  my $stats = {};
  for my $line (split /\n/, $stdout) {
    my ($attr, $value) = split(/\t/, $line);
    $attr =~ s/\:$//;
    if ($attr eq "raw total sequences" 
      || $attr eq "reads properly paired" 
      || $attr eq "reads mapped" 
      || $attr eq "average length") {
      $stats->{$attr} = $value;
    }
  }

  return $stats;
}

1;
