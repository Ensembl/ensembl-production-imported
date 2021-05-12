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

package Bio::EnsEMBL::EGPipeline::DNAFeatures::RepeatMasker;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable::RepeatMasker;
use File::Path qw(make_path);

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::AnalysisRun');

sub fetch_runnable {
  my ($self) = @_;
  
  my %parameters;
  if (%{$self->param('parameters_hash')}) {
    %parameters = %{$self->param('parameters_hash')};
  }
  
  #my $repeatmasker_cache = $self->param('repeatmasker_cache');
  #if (defined $repeatmasker_cache && -e $repeatmasker_cache) {
  #  $ENV{'REPEATMASKER_CACHE'} = $repeatmasker_cache;
  #}
  
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::RepeatMasker->new
  (
    -query    => $self->param('query'),
    -program  => $self->param('program'),
    -analysis => $self->param('analysis'),
    -datadir  => $self->param('datadir'),
    -bindir   => $self->param('bindir'),
    -libdir   => $self->param('libdir'),
    -workdir  => $self->param('workdir'),
    -timer    => $self->param('timer'),
    %parameters,
  );
  
  $self->param('save_object_type', 'RepeatFeature');
  
  return $runnable;
}

sub results_by_index {
  my ($self, $results) = @_;
  my %seqnames;
  
  my ($header, $body) = $results =~ /(.+\n\n)(.+)/ms;
  my @lines = split(/\n/, $body);
  foreach my $line (@lines) {
    my ($seqname) = $line =~ /^\s*(?:\S*\s+){4}(\S+)/;
    $seqnames{$seqname}{'result'} .= "$line\n";
  }
  foreach my $seqname (keys %seqnames) {
    $seqnames{$seqname}{'header'} = $header;
  }
  
  return %seqnames;
}

1;
