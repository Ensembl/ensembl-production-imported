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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::InferStrandness;

use strict;
use warnings;
use Capture::Tiny ':all';
use Try::Tiny;
use File::Spec::Functions qw(catdir);

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');


sub param_defaults {
  my ($self) = @_;
  
  return {
    keep_ambiguous => 1,  # If ambiguous, keep with a flag
    use_input_if_ambiguous => 0, # Or, use the input if any
    n_samples => 200000,
    infer_max => 0.80,  # Above this, infer stranded
    infer_min => 0.65,  # Between this and infer_max, ambiguous
    infer_failed_max => 0.2, # Threshold for max failed reads
    infer_mean_reads => 1_000, # Minimum number of usable reads to make an inference
  };
}

sub run {
  my ($self) = @_;
  
  my $bed_file = $self->param_required('bed_file');
  my $sam_file = $self->param_required('sam_file');
  my $n_samples = $self->param_required('n_samples');
  
  my $input_metadata = $self->param('input_metadata');
  my $input_is_paired = $input_metadata->{'is_paired'};
  my $input_is_stranded = $input_metadata->{'is_stranded'};

  my $sample_1 = $self->param('sample_seq_file_1');
  my $sample_2 = $self->param('sample_seq_file_2');

  my $res_dir = $self->param_required('results_dir');
  my $log_file = catdir($res_dir, 'log.txt');
  
  my ($strandness, $is_paired, $is_ambiguous, $no_reads, $output) = $self->get_strandness($bed_file, $sam_file, $n_samples, $input_is_stranded);

  my $is_stranded = 0;
  my $strand_direction = "";
  if ($strandness) {
    $is_stranded = 1;

    if ($strandness eq 'RF' or $strandness eq 'R') {
      $strand_direction = 'reverse';
    } elsif ($strandness eq 'FR' or $strandness eq 'F') {
      $strand_direction = 'forward';
    }
  }

  open my $LOG, ">>", $log_file;

  print $LOG "Strandness inference for '".$self->param('sample_name')."': $output\n";

  # Check inferred vs input
  if (defined $input_is_paired and $input_is_paired != $is_paired) {
    my $input = "Input = " . ($input_is_paired ? "paired" : "single");
    my $infer = "Infer = " . ($is_paired ? "paired" : "single");
    print $LOG "WARNING: input and inferred paired-end differ: $input vs $infer\n";
  }
  if (defined $input_is_stranded and $input_is_stranded != $is_stranded) {
    my $input = "Input = " . ($input_is_stranded ? "strand-specific" : "unstranded");
    my $infer = "Infer = " . ($is_stranded ? "strand-specific" : "unstranded");
    print $LOG "WARNING: input and inferred strand-specificity differ: $input vs $infer\n";
  }
  if (defined $strand_direction) {
    print $LOG "Inferred direction = $strand_direction\n";
  }
  if (not defined $is_paired) {
    if (defined $input_is_paired) {
      $is_paired = $input_is_paired;
      print $LOG "WARNING: could not infer paired/single-end status. Using input value: $is_paired\n";
    } elsif ($no_reads) {
      $is_paired = ($sample_1 and $sample_2) ? 1 : 0;
      print $LOG "WARNING: could not infer paired/single-end status: no useful reads in the sample. Using value based on number of files: $is_paired\n";
    } else {
      print $LOG "WARNING: could not infer paired/single-end status\n";
      die("Could not infer paired/single-end status");
    }
  }
  if (not defined $is_stranded) {
    if (defined $input_is_stranded) {
      print $LOG "WARNING: could not infer strand-specificity. Using input value: $input_is_stranded\n";
      $is_stranded = $input_is_stranded;
    } elsif ($no_reads) {
      $is_stranded = 0;
      print $LOG "WARNING: could not infer strand-specificity. Using default value: $is_stranded\n";
    } else {
      print $LOG "WARNING: could not infer strand-specificity\n";
      die("Could not infer strand-specificity");
    }
  }
  if ($is_ambiguous) {
    print $LOG "WARNING: This inference in ambiguous\n";
  }
  close $LOG;

  my $aligner_metadata = {
      is_stranded => $is_stranded,
      is_paired => $is_paired,
      strand_direction => $strand_direction,
  };
  if ($is_ambiguous) {
    $aligner_metadata->{stranded_ambiguous} = 1;
  }
  $self->dataflow_output_id({
      aligner_metadata => $aligner_metadata
    },  2);

  cleanup_file($sam_file);
  cleanup_file($self->param('sample_seq_file_1'));
  cleanup_file($self->param('sample_seq_file_2'));
}

sub get_strandness {
  my ($self, $bed_file, $sam_file, $n_samples, $input_is_stranded) = @_;

  # Default empty = unstranded
  my $strandness = '';

  # Prepare the command for the RSeQ script
  my $cmd = "infer_experiment.py";
  my @cmd_args = (
    "-r $bed_file",
    "-i $sam_file",
    "-s $n_samples",
  );
  $cmd .= " " . join( " ", @cmd_args );
  
  try {
    my ($stdout, $stderr, $exit) = capture {
      system( $cmd );
    };

    die("No inference output: $stderr") if not $stdout;
    
    # Check number of usable reads
    my $n_reads = 0;
    if ($stderr =~ /Total (\d+) usable reads were sampled/) {
      $n_reads = $1;
    }
    my $min_reads = $self->param('infer_mean_reads');
    if ($n_reads < $min_reads) {
      die("Inference failed (number of usable reads is too small: $n_reads < $min_reads): $stderr");
    }
    
    # Check success
    if ($exit != 0) {
      die("Inference failed ($exit): $stderr");
    } else {
      my $output = "$stdout\n$stderr";
      my ($strandness, $is_paired, $is_ambiguous) = $self->parse_inference($stdout, $stderr, $input_is_stranded);
      return ($strandness, $is_paired, $is_ambiguous, $n_reads, $output);
    }
  } catch {
    # Nothing could be determined? Return nothing but continue
    die("Inference failed with command:\n$cmd\nOutput: $_");
  };
}

sub cleanup_file {
  my ($file) = @_;

  unlink $file if $file;
}

sub parse_inference {
  my ($self, $text, $err, $input_is_stranded) = @_;

  my $is_paired;
  my %strand = (
    '++,--' => 'forward',
    '+-,-+' => 'reverse',
    '1++,1--,2+-,2-+' => 'forward',
    '1+-,1-+,2++,2--' => 'reverse',
  );

  my %stats;
  my $is_ambiguous = 0;

  for my $line (split /[\r\n]+/, $text) {
    # Single or paired-end
    if ($line =~ /This is (.+) Data/) {
      my $lib = $1;
      if ($lib eq 'SingleEnd') {
        $is_paired = 0;
      } elsif ($lib eq 'PairEnd') {
        $is_paired = 1;
      }
    } elsif ($line =~ /^Fraction of reads failed to determine: (.+)$/) {
      $stats{failed} = $1;
    } elsif ($line =~ /^Fraction of reads explained by "(.+?)": (.+)$/) {
      my ($st, $fraction) = ($1, $2);

      if (exists $strand{$st}) {
        $stats{ $strand{$st} } = $fraction + 0;
      }
    }
    elsif ($line =~ /Unknown Data type/) {
      die("Unknown data type!: (stderr: $err)");
    }
  }
  $stats{failed} //= 0;

  if (not defined $is_paired) {
    die("Could not determine the library: single or paired-end from:\n$text");
  }

  if (not $stats{ forward }) {
    warn("Could not parse forward data");
  }
  if (not $stats{ reverse }) {
    warn("Could not parse reverse data");
  }

  # Default: unstranded
  my $strandness = '';

  my $aligned = (1 - $stats{failed});
  my $max_failed = $self->param('infer_failed_max');
  my $max_ambiguous = $aligned * $self->param('infer_max'); # Anything above is considered stranded
  my $min_ambiguous = $aligned * $self->param('infer_min');
  
  # Too much failed: can't infer strandness
  if ($stats{failed} > $max_failed) {
    $strandness = '';
    die("Run contains too many reads failed to align ($stats{failed}). Please check: $text");

  # Stranded forward
  } elsif ($stats{ forward } > $min_ambiguous) {
    $strandness = $is_paired ? 'FR' : 'F';
    
    # Not enough power to infer: die or use input
    if ($stats{ forward } < $max_ambiguous and not $input_is_stranded) {
      if ($self->param('keep_ambiguous')) {
        $is_ambiguous = 1;
      } elsif ($self->param('use_input_if_ambiguous')) {
        $strandness = '';
      } else {
        die("Ambiguous strand inference: $stats{ forward } < $max_ambiguous:\n$text");
      }
    }

  # Stranded reverse
  } elsif ($stats{ reverse } > $min_ambiguous) {
    $strandness = $is_paired ? 'RF' : 'R';

    # Not enough power to infer: use input
    if ($stats{ reverse } < $max_ambiguous and not $input_is_stranded) {
      if ($self->param('keep_ambiguous')) {
        $is_ambiguous = 1;
      } elsif ($self->param('use_input_if_ambiguous')) {
        $strandness = '';
      } else {
        die("Ambiguous strand inference: $stats{ reverse } < $max_ambiguous:\n$text");
      }
    }
  }

  return ($strandness, $is_paired, $is_ambiguous);
}

1;
