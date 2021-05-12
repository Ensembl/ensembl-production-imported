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

package Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::EmailReport;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EmailReport');

use Bio::EnsEMBL::EGPipeline::Common::Aligner;

sub param_defaults {
  my ($self) = @_;
  
  return {
    'samtools_dir' => undef,
  };
}

sub fetch_input {
	my ($self) = @_;
  
  my $species      = $self->param_required('species');
  my $samtools_dir = $self->param('samtools_dir');
  my $bam_file     = $self->param_required('merged_bam_file');
  my $bw_file      = $self->param('bw_file');
  my $ini_file     = $self->param_required('ini_file');
  my $cmds_file    = $self->param_required('cmds_file');
  
  my $aligner = Bio::EnsEMBL::EGPipeline::Common::Aligner->new(
    -samtools_dir => $samtools_dir,
  );
  my $stats = $aligner->get_bam_stats($bam_file);
  
  my $text = 
    "The DNA Sequence Alignment pipeline for $species has generated ".
    "a BAM file ($bam_file) and index, ";
  
  if (defined $bw_file && -e $bw_file) {
    $text .= "a BigWig file ($bw_file), ";
  }
  
  $text .= 
    "and an accompanying INI file ($ini_file).\n\n".
    "The commands for generating these results are ".
    "stored in a JSON file ($cmds_file)\n\n".
    "BAM statistics:\n$stats\n\n";
  
  $self->param('text', $text);
}

1;
