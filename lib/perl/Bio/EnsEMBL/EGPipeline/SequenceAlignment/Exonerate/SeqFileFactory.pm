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

package Bio::EnsEMBL::EGPipeline::SequenceAlignment::Exonerate::SeqFileFactory;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Path::Tiny qw(path);

sub write_output {
  my ($self) = @_;
  
  my $species          = $self->param_required('species');
  my $seq_file         = $self->param('seq_file');
  my $seq_file_species = $self->param('seq_file_species');
  my $data_type        = $self->param('data_type');
  my $reformat_header  = $self->param('reformat_header');
  my $trim_est         = $self->param('trim_est');
  
  my $flow = $self->param_required('use_exonerate_server') ? 3 : 2;
  
  if (exists $$seq_file_species{$species}) {
    push @$seq_file, $$seq_file_species{$species};
  }
  
  my $counter = 1;
  
  foreach my $fasta_file (@$seq_file) {
    if ($data_type eq 'est' && $trim_est) {
      $fasta_file = $self->trim_est($fasta_file);
    }
    
    if ($reformat_header) {
      $self->reformat_header($fasta_file);
    }
    
    my $dataflow_output = {
      'fasta_file'  => $fasta_file,
      'out_subdir' => 'dataset_'.$counter++,
    };
    $self->dataflow_output_id($dataflow_output, $flow);
  }
}

sub reformat_header {
  my ($self, $fasta_file) = @_;
  
  my $file = path($fasta_file);
  my $data = $file->slurp;
  $data =~ s/^>gi\|\d+\|gb\|([^\|]+)\|\S+/>$1/g;
  $file->spew($data);
}

sub trim_est {
  my ($self, $fasta_file) = @_;
  
  my $trimest_exe = $self->param('trimest_exe');
  (my $trimmed_file = $fasta_file) =~ s/(\.\w+)$/-trimmed$1/;
  
  my $cmd = "$trimest_exe -seq $fasta_file -out $trimmed_file";
  system($cmd) == 0 || $self->throw("Cannot execute $cmd");
  
  return $trimmed_file;
}

1;
