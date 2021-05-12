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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::ReadMetadata;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Capture::Tiny ':all';
use Path::Tiny qw(path);
use File::Basename qw(dirname);
use File::Spec::Functions qw(catdir);
use JSON;

sub run {
  my ($self) = @_;

  my $results_dir = $self->param_required('results_dir');
  my $meta_file = catdir($results_dir, 'metadata.json');

  # Slurp the json file
  open my $meta_fh, '<', $meta_file;
  $/ = undef;
  my $data_json = <$meta_fh>;
  close $meta_fh;
  my $data = decode_json($data_json);
  
  # Get the data that we need for HTSeq
  my %metadata = ();
  my %map_bool = (
    "is_paired" => "hasPairedEnds",
    "is_stranded" => "isStrandSpecific",
  );
  my %map_value = (
    "strand_direction" => "strandDirection",
  );
  
  for my $key (keys %map_value) {
    my $value = $map_value{$key};
    $metadata{$key} = $data->{$value} if exists $data->{$value};
  }
  for my $key (keys %map_bool) {
    my $value = $map_bool{$key};
    $metadata{$key} = $data->{$value} ? 1 : 0 if exists $data->{$value};
  }

  $self->dataflow_output_id(\%metadata, 2);
}

1;

