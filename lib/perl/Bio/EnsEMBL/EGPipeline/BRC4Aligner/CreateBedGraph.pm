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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::CreateBedGraph;

use strict;
use warnings;
use Capture::Tiny ':all';
use base ('Bio::EnsEMBL::EGPipeline::BRC4Aligner::Base');

sub param_defaults {
  my ($self) = @_;
  
  return {
    'bedtools_dir'  => undef,
    'ucscutils_dir' => undef,
    'clean_up'      => 1,
  };
}

sub run {
	my ($self) = @_;
  
  my $bamutils_dir  = $self->param('bamutils_dir');
  my $bam_file      = $self->param_required('bam_file');
  my $strand  = $self->param('strand');
  my $direction  = $self->param('strand_direction');
  
  my $bamutils = 'bamutils';
  if (defined $bamutils_dir) {
    $bamutils = "$bamutils_dir/$bamutils";
  }

  my ($bed_file, $bed_cmd) = $self->convert_to_bed($bam_file, $strand, $direction);
  
  $self->param('bed_file', $bed_file);
  $self->param('cmds',    $bed_cmd);
}

sub write_output {
  my ($self) = @_;

    my $align_cmds = {
      cmds => $self->param('cmds'),
      sample_name => $self->param('sample_name'),
    };
    $self->store_align_cmds($align_cmds);
  
  my $dataflow_output = {
    'bed_file' => $self->param('bed_file'),
  };
  
  $self->dataflow_output_id($dataflow_output, 2);
}

sub convert_to_bed {
  my $self = shift;
  my ($bam, $strand, $direction) = @_;

  my $bed = $bam;
  $bed =~ s/\.bam$/.bed/;

  my $filter = "";
  
  # Even a strand is given, make sure we are on the right one
  if ($strand and $direction) {
    if ($strand eq 'firststrand') {
      if ($direction eq 'forward') {
        $filter = "-plus";
      } else {
        $filter = "-minus";
      }
    } elsif ($strand eq 'secondstrand') {
      if ($direction eq 'forward') {
        $filter = "-minus";
      } else {
        $filter = "-plus";
      }
    }
  } else {
    print("No strand used.\n");
  }

  my $bed_tmp = "$bed.unsorted";
  my $cmd = "bamutils tobedgraph $filter $bam > $bed_tmp && bedSort $bed_tmp $bed";
  
  my ($stdout, $stderr, $exit) = capture {
    system($cmd);
  };
  if ($exit) {
    die("Cannot execute $cmd:\n$stderr");
  }
  unlink $bed_tmp;
  
  # Check the bed file is not empty
  if (not -s $bed) {
    $self->throw("Bed file '$bed' is empty. Made with strand='$strand' and direction='$filter'. Command: $cmd");
  }

  return ($bed, $cmd);
}

1;
