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

Bio::EnsEMBL::EGPipeline::ProjectGeneDesc::SourceFactory

=head1 DESCRIPTION

Parse the config hash that defines what to project to what.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::ProjectGeneDesc::SourceFactory;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');
use File::Spec::Functions qw(catdir);

sub write_output {
  my ($self)  = @_;
  
  my $config     = $self->param_required('config');
  my $flow       = $self->param_required('flow');
  my $output_dir = $self->param_required('output_dir');
  
  foreach my $id (sort { $$flow{$a} <=> $$flow{$b} } keys %$flow) {
    my $output_ids = $$config{$id};
    $$output_ids{'projection_dir'} = catdir($output_dir, $id);
    
    $self->dataflow_output_id($output_ids, $$flow{$id});
  }
}

1;
