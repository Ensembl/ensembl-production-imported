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

package Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::DumpGTF;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::FileDump::GTFDumper');

sub run {
  my ($self) = @_;
  
  my $out_file = $self->param_required('out_file');
  
  if (!-e $out_file) {
    $self->SUPER::run;
  }
}

1;
