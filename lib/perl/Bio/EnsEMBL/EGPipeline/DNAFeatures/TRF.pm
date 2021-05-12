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

package Bio::EnsEMBL::EGPipeline::DNAFeatures::TRF;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable::TRF;

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::AnalysisRun');

sub fetch_runnable {
  my ($self) = @_;
  
  my %parameters;
  if (%{$self->param('parameters_hash')}) {
    %parameters = %{$self->param('parameters_hash')};
  }
  
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::TRF->new
  (
    -query    => $self->param('query'),
    -program  => $self->param('program'),
    -analysis => $self->param('analysis'),
    -datadir  => $self->param('datadir'),
    -bindir   => $self->param('bindir'),
    -libdir   => $self->param('libdir'),
    -workdir  => $self->param('workdir'),
    %parameters,
  );
  
  $self->param('save_object_type', 'RepeatFeature');
  
  return $runnable;
}

sub update_options {
  my ($self, $runnable) = @_;
  
  # The 'options' have further options appended when the command is run.
  my $options = $runnable->options." -d -h";
  
  my $analysis = $self->param('analysis');
  if ($analysis->parameters ne $options) {
    $analysis->parameters($options);
    $self->param('analysis_adaptor')->update($analysis);
  }
}

sub results_by_index {
  my ($self, $results) = @_;
  my %seqnames;
  
  my @results = split(/Sequence:/, $results);
  my $header = shift @results;
  foreach my $result (@results) {
    my ($seqname) = $result =~ /^\s*(\S+)/;
    $seqnames{$seqname}{'result'} = "Sequence:$result";
    $seqnames{$seqname}{'header'} = $header;
  }
  
  return %seqnames;
}

sub filter_output {
  my ($self, $runnable) = @_;
  my %deduplicated;
  
  # There are sometimes duplicate rows, because only certain columns are
  # taken from the raw output; so if the only difference between two results
  # is in the columns that are _not_ used, when the data is extracted those
  # rows end up being duplicates.
  foreach my $rf (@{$runnable->output}) {
    my $key = join(",", ($rf->start, $rf->end, $rf->repeat_consensus->seq, $rf->score));
    if (! exists $deduplicated{$key}) {
      $deduplicated{$key} = $rf;
    } else {
      $self->warning('TRF duplicate on sequence '.$runnable->query->name." for key $key.");
    }
  }
      
  my @deduplicated = values %deduplicated;
  # Output is cumulative, so need to manually erase existing results.
  # (Note that calling the runnable's 'output' method will NOT work.
  $runnable->{'output'} = \@deduplicated;
}

1;
