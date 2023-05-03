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
    force_aligner_metadata => undef,
    ignore_single_paired => 0,
    ambiguity_majority_rule => 0.10, # Ratio of acceptable ambiguous samples
  };
}

sub fetch_input {
  my ($self) = @_;
}

sub run {
  my ($self) = @_;
  
  my $force_metadata = $self->param('force_aligner_metadata');
  
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
  if (not $self->param('ignore_single_paired') and not defined $met->{is_paired}) {
    push @errors, "No 'is_paired' value";
  }
  push @errors, "No 'is_stranded' value" if not defined $met->{is_stranded};
  if ($met->{is_stranded}) {
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
  );
  
  my @stranded_ambiguous;
  for my $sample (keys %$metadata_hash) {
    my $met = $metadata_hash->{$sample};
    for my $key (keys %$met) {
      my $value = $met->{$key};
      if ($key eq 'stranded_ambiguous' and $value == 1) {
        push @stranded_ambiguous, $met;
      } else {
        $met_counts{$key}{$value}++;
      }
    }
  }
  
  # Special: we can process mixed single/paired end
  if ($self->param('ignore_single_paired')) {
    delete $met_counts{'is_paired'};
  }

  # Check if there are ambiguous samples to check
  my $ambiguity_ratio = 0;
  if ($self->param('ambiguity_majority_rule') and @stranded_ambiguous) {
    my $num_ambiguous = scalar @stranded_ambiguous;
    my $num_samples = scalar keys %$metadata_hash;
    $ambiguity_ratio = $num_ambiguous / $num_samples;
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

  # Compare ambiguity
  if ($self->param('ambiguity_majority_rule') and $ambiguity_ratio) {
    if ($ambiguity_ratio > $self->param('ambiguity_majority_rule')) {
      push @errors, "Too many ambiguous samples: $ambiguity_ratio";
    } else {
      # Make sure that the ambiguous samples are compatible with the consensus
      for my $ambig (@stranded_ambiguous) {
        if ($ambig->{strand_direction} != $consensus{strand_direction}) {
          push @errors, "Ambiguous sample strand_direction differs from consensus: " . $ambig->{strand_direction} . " != " . $consensus{strand_direction};
        }
      }
    }
  }
  
  if (@errors) {
    my $json = JSON->new->pretty->canonical(1);
    my $cons_str = $json->encode(\%consensus);
    $cons_str =~ s/:/ => /g;
    
    my $all_met_str = $json->encode($metadata_hash);
    die("Could not create a consensus: " . join("; ", @errors) . "\nPlease check the logs and fix the Current Consensus (and copy it as 'force_aligner_metadata' param for this job):\n" . $cons_str . "\nAll metadata: " . $all_met_str);
  }

  return \%consensus;
}

1;
