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

Bio::EnsEMBL::EGPipeline::Common::ParameterFromJson

=head1 DESCRIPTION

This simple module loads parameters from a json file. The json file is expected
to be a dictionary for which the keys and their values will become the
parameters names and values.

If a json file is given, the parameters will be output to dataflow 1.
So whether a file is given or not, the pipeline will continue.

=head1 Author

Matthieu Barba

=cut

package Bio::EnsEMBL::EGPipeline::Common::RunnableDB::ParametersFromJson;

use strict;
use warnings;
use JSON;

use base (
  'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base'
);

sub run {
  my $self = shift @_;
  
  my $json_path = $self->param('parameters_json');
  
  if ($json_path) {
    my $data = $self->decode_json_file($json_path);
    $self->dataflow_output_id($data, 1);
  }
}

sub decode_json_file {
  my ($self, $path) = @_;

  my $json;
  {
    local $/; #enable slurp
    open my $fh, "<", $path;
    $json = <$fh>;
  }

  my $data = decode_json($json);
  return $data;
}

1;
