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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::AggregateMetadata;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use JSON;

sub param_defaults {
  my ($self) = @_;
  
  return {
  };
}

sub fetch_input {
  my ($self) = @_;
}

sub run {
  my ($self) = @_;
  
  my $all_metadata = $self->param('aligner_metadata_array');
  my $consensus_metadata = $self->create_consensus_metadata($all_metadata);

  my $output_data = {
    aggregated_aligner_metadata => $consensus_metadata,
  };
  $self->dataflow_output_id($output_data, 2);
}

sub create_consensus_metadata {
  my ($self, $metadata_array) = @_;

  my %met_counts = (
    "is_paired" => {},
    "is_stranded" => {},
    "strand_direction" => {},
    "strandness" => {}
  );
  
  for my $met (@$metadata_array) {
    for my $key (keys %$met) {
      my $value = $met->{$key};
      $met_counts{$key}{$value}++;
    }
  }
  
  # Check that there is only 1 value for each key
  my %consensus;
  my @errors;
  for my $key (keys %met_counts) {
    my @sub_keys = keys %{$met_counts{$key}};
    
    if (@sub_keys == 1) {
      my $unique_value = shift @sub_keys;
      $consensus{$key} = $unique_value;
    } else {
      push @errors, "There are " . scalar(@sub_keys) . " different values for metadata $key.";
    }
  }
  
  if (@errors) {
    die("Could not create a consensus: " . join("; ", @errors));
  }

  return \%consensus;
}

1;
