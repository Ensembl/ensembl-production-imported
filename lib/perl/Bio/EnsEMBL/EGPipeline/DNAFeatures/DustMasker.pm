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

package Bio::EnsEMBL::EGPipeline::DNAFeatures::DustMasker;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable::DustMasker;

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::AnalysisRun');

sub fetch_runnable {
  my $self = shift @_;
  
  my %parameters;
  if (%{$self->param('parameters_hash')}) {
    %parameters = %{$self->param('parameters_hash')};
  }
  
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::DustMasker->new
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
  
  # The 'options' are updated when the command is run.
  my $options = $runnable->options;
  
  my $analysis = $self->param('analysis');
  if ($analysis->parameters ne $options) {
    $analysis->parameters($options);
    $self->param('analysis_adaptor')->update($analysis);
  }
}

sub results_by_index {
  my ($self, $results) = @_;
  my %seqnames;
  
  my @results = split(/>/, $results);
  my $header = shift @results;
  foreach my $result (@results) {
    my ($seqname) = $result =~ /^\s*(\S+)/;
    $seqnames{$seqname}{'result'} = ">$result";
    $seqnames{$seqname}{'header'} = $header;
  }
  
  return %seqnames;
}

1;
