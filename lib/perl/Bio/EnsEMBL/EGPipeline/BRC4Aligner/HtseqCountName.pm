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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::HtseqCountName;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::BRC4Aligner::Base');

use Capture::Tiny ':all';
use Path::Tiny qw(path);
use File::Basename qw(dirname);
use File::Spec::Functions qw(catfile);

sub param_defaults {
  my ($self) = @_;
  
  return {
    # Features on which to count the read
    features => 'exon',
  };
}

sub run {
  my ($self) = @_;

  $self->dbc and $self->dbc->disconnect_if_idle();

  my $bam = $self->param_required('bam_file');
  my $gtf = $self->param_required('gtf_file');
  my $strand = $self->param_required('strand');
  my $number = $self->param_required('number');
  my $strand_direction = $self->param('strand_direction');
  my $feature = $self->param('feature');

  my $results_dir = dirname($bam);

  my $htseq_file = 'genes.htseq-union.' . $feature . '.' . $strand;
  if ($number eq 'total') {
    $htseq_file .= '.nonunique';
  }
  #$htseq_file .= '.' . $feature;
  $htseq_file .= '.counts';

  $htseq_file = catfile($results_dir, $htseq_file);
  
  print("Creating htseqcount file $htseq_file");

  # Define stranded and nonunique flags
  my $number_flag = "";
  if ($number eq 'total') {
    $number_flag = "--nonunique all";
  } elsif ($number eq 'unique') {
    $number_flag = "--nonunique none";
  } else {
    die("Unrecognized number: $number");
  }
  my $stranded_flag = "";

  if ($strand ne 'unstranded') {
    die "No strand direction given for stranded data" if not $strand_direction;
    
    if ($strand eq 'firststrand' and $strand_direction eq 'forward'
        or $strand eq 'secondstrand' and $strand_direction eq 'reverse') {
      $stranded_flag = "--stranded=yes";
    } elsif ($strand eq 'firststrand' and $strand_direction eq 'reverse'
        or $strand eq 'secondstrand' and $strand_direction eq 'forward') {
      $stranded_flag = "--stranded=reverse";
    } else {
      die("Unrecognized strand: $strand, $strand_direction");
    }
  } else {
    $stranded_flag = "--stranded=no";
  }
  my $feature_flag;
  if ($feature) {
    $feature_flag = "--type=$feature";
  }
  my $params = "$number_flag $stranded_flag $feature_flag";

  # Run
  my $cmd = $self->run_htseq_count($bam, $gtf, $htseq_file, $params);
  
  $self->param("cmd", $cmd);
}

sub write_output {
  my ($self) = @_;
  
  my $cmd = $self->param("cmd");
  
  my $version = $self->get_htseq_version();

  my $align_cmds = {
    cmds => $cmd,
    sample_name => $self->param('sample_name'),
    'version'  => $version,
  };
  $self->store_align_cmds($align_cmds);
  
  $self->dataflow_output_id({ cmds => $cmd },  2);
}

sub get_htseq_version {
  my ($self) = @_;

  my $cmd = "htseq-count --version a b";
  my ($stdout, $stderr, $exit) = capture {
    system($cmd);
  };

  if ($exit != 0) {
    die("Could not get htseq-count version: $stderr");
  }
  
  my $version = chomp $stdout;

  return $version;
}

sub run_htseq_count {
  my ($self, $bam, $gtf, $htseq_file, $params, $by_pos) = @_;

  my $order = $by_pos ? '--order=pos' : '--order=name';

  my $base_cmd = "htseq-count -a 0 --format=bam --idattr=gene_id --mode=union -q";
  my $cmd = "$base_cmd $order $params $bam $gtf > $htseq_file";
  my ($stdout, $stderr, $exit) = capture {
    system($cmd);
  };

  if ($exit != 0) {
    if ($stderr =~ /'NoneType' object has no attribute 'encode/) {
      die("Index failed (special case): $stderr");
    } else {
      die("Index failed: $stderr");
    }
  }

  return $cmd;
}

1;
