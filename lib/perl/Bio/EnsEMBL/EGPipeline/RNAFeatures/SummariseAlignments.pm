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


=pod

=head1 NAME

Bio::EnsEMBL::EGPipeline::RNAFeatures::SummariseAlignments

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::RNAFeatures::SummariseAlignments;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub run {
  my ($self) = @_;
  my $run_cmscan   = $self->param_required('run_cmscan');
  my $run_trnascan = $self->param_required('run_trnascan');
  my $load_mirbase = $self->param_required('load_mirbase');
  my $pipeline_dir = $self->param_required('pipeline_dir');
  
  if ($run_cmscan) {
    my $cmd = "cat $pipeline_dir/*/cmscan.txt > $pipeline_dir/cmscan.txt";
    system($cmd) == 0 || $self->throw("Failed to execute $cmd");
  }
  
  if ($run_trnascan) {    
    my $cmd = "cat $pipeline_dir/*/trnascan.txt > $pipeline_dir/trnascan.txt";
    system($cmd) == 0 || $self->throw("Failed to execute $cmd");
  }
  
  if ($load_mirbase) {
    my $cmd = "cat $pipeline_dir/*/mirbase.txt > $pipeline_dir/mirbase.txt";
    system($cmd) == 0 || $self->throw("Failed to execute $cmd");
  }
  
}

sub write_output {
  my ($self) = @_;
  my $run_cmscan   = $self->param_required('run_cmscan');
  my $pipeline_dir  = $self->param_required('pipeline_dir');
  my $evalue_levels = $self->param_required('evalue_levels');
  
  if ($run_cmscan) {
    foreach my $evalue (keys %$evalue_levels) {
      my $output_ids = {
        'cmscanfile'   => "$pipeline_dir/cmscan.txt",
        'evalue'       => $evalue,
        'biotypesfile' => "$pipeline_dir/biotypes_$evalue.svg",
        'distinctfile' => "$pipeline_dir/distinct_$evalue.svg",
        'plotcolour'   => $$evalue_levels{$evalue},
      };
      
      $self->dataflow_output_id($output_ids, 2);
    }
  }
}

1;
