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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::CreateBigWig;

use strict;
use warnings;
use Capture::Tiny ':all';
use base ('Bio::EnsEMBL::EGPipeline::BRC4Aligner::Base');

sub param_defaults {
  my ($self) = @_;
  
  return {
    'bedtools_dir'  => undef,
    'ucscutils_dir' => undef,
  };
}

sub run {
	my ($self) = @_;
  
  my $bedtools_dir  = $self->param('bedtools_dir');
  my $ucscutils_dir = $self->param('ucscutils_dir');
  my $length_file   = $self->param_required('length_file');
  my $bam_file      = $self->param_required('bam_file');
  
  my $wig_file = "$bam_file.wig";
  (my $bw_file = $bam_file) =~ s/\.bam$/\.bw/;
  if ($bw_file eq $bam_file) {
    $bw_file = "$bam_file.bw";
  }
  
  my $bedtools = 'bedtools';
  if (defined $bedtools_dir) {
    $bedtools = "$bedtools_dir/$bedtools"
  }
  my $ucscutils = 'wigToBigWig';
  if (defined $ucscutils_dir) {
    $ucscutils = "$ucscutils_dir/$ucscutils"
  }
  
  my $wig_cmd =
    "$bedtools genomecov ".
    " -g $length_file ".
    " -ibam $bam_file ".
    " -bg ".
    " -split ".
    " > $wig_file ";
  my $bw_cmd =
    "$ucscutils ".
    " $wig_file ".
    " $length_file ".
    " $bw_file ";
  
  # Reuse precalculated bigwig if it was finished
  if (not -s $bw_file or -s $wig_file) {
    $self->_execute($wig_cmd);
    $self->_execute($bw_cmd);

    unlink $wig_file;
  }
  
  $self->param('wig_cmd', $wig_cmd);
  $self->param('bw_cmd', $bw_cmd);
}

sub write_output {
  my ($self) = @_;
  
  $self->store_align_cmds({
    cmds => $self->param('wig_cmd'),
    sample_name => $self->param('sample_name'),
  });
  $self->store_align_cmds({
    cmds => $self->param('bw_cmd'),
    sample_name => $self->param('sample_name'),
  });
}

sub _execute {
  my $self = shift;
  my ($cmd) = @_;
  
  
  my ($stdout, $stderr, $exit) = capture {
    system($cmd);
  };
  if ($exit) {
    $self->throw("Cannot execute $cmd:\n$stderr");
  }
}

1;
