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
    force_aligner_metadata => undef
  };
}

sub fetch_input {
  my ($self) = @_;
}

sub run {
  my ($self) = @_;
  
  
  my $force_metadata = $self->param('force_metadata');
  
  my $consensus_metadata;
  if ($force_metadata) {
    $consensus_metadata = $force_metadata;
  } else {
    my $all_metadata = $self->param('aligner_metadata_hash');
    $consensus_metadata = $self->create_consensus_metadata($all_metadata);
  }

  # Check the metadata
  $self->check_metadata($consensus_metadata);

  my $output_data = {
    aggregated_aligner_metadata => $consensus_metadata,
  };
  $self->dataflow_output_id($output_data, 2);
}

sub check_metadata {
  my ($self, $met) = @_;
  
  my @errors;
  push @errors, "No 'is_paired' value" if not defined $met->{is_paired};
  push @errors, "No 'is_stranded' value" if not defined $met->{is_stranded};
  if ($met->{is_stranded}) {
    push @errors, "No 'strandness' value" if not defined $met->{strandness};
    push @errors, "No 'strand_direction' value" if not defined $met->{strand_direction};
  }
  
  die("Error in consensus metadata: " . join("\n", @errors)) if @errors;
}
  
sub create_consensus_metadata {
  my ($self, $metadata_hash) = @_;

  my %met_counts = (
    "is_paired" => {},
    "is_stranded" => {},
    "strand_direction" => {},
    "strandness" => {}
  );
  
  for my $sample (keys %$metadata_hash) {
    my $met = $metadata_hash->{$sample};
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
    @sub_keys = grep { $_ ne "" } @sub_keys;
    
    if (@sub_keys > 1) {
      my @pairs = map { "$_=$met_counts{$key}->{$_}" } @sub_keys;
      push @errors, "There are " . scalar(@sub_keys) . " different values for metadata '$key' (".join(", ", @pairs).")";
      $consensus{$key} = join("|", @sub_keys);
    } else {
      my $unique_value = shift @sub_keys;
      $consensus{$key} = $unique_value;
    }
  }
  
  if (@errors) {
    my $cons_str = encode_json(\%consensus);
    $cons_str =~ s/:/ => /g;
    
    my $all_met_str = encode_json($metadata_hash);
    die("Could not create a consensus: " . join("; ", @errors) . "\nPlease check the logs and fix the Current Consensus (and copy it as 'force_metadata' param for this job):\n" . $cons_str . "\nAll metadata: " . $all_met_str);
  }

  return \%consensus;
}

1;
