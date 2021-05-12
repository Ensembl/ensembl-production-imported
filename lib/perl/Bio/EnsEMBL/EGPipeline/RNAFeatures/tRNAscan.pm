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

package Bio::EnsEMBL::EGPipeline::RNAFeatures::tRNAscan;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable::tRNAscan;
use File::Path qw(make_path);

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::AnalysisRun');

sub param_defaults {
  my $self = shift @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'pseudo'    => 0,
    'threshold' => 20,
  };
}

sub fetch_runnable {
  my ($self) = @_;
  
  my %parameters;
  if (%{$self->param('parameters_hash')}) {
    %parameters = %{$self->param('parameters_hash')};
  }
  
  my $db_name = $self->param_required('db_name');
  my $external_db_id = $self->fetch_external_db_id($db_name);
  $self->throw("Unrecognised external DB '$db_name'") unless $external_db_id;
  
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::tRNAscan->new
  (
    -query          => $self->param('query'),
    -program        => $self->param('program'),
    -analysis       => $self->param('analysis'),
    -datadir        => $self->param('datadir'),
    -bindir         => $self->param('bindir'),
    -libdir         => $self->param('libdir'),
    -workdir        => $self->param('workdir'),
    -external_db_id => $external_db_id,
    -pseudo         => $self->param('pseudo'),
    -threshold      => $self->param('threshold'),
    %parameters,
  );
  
  $self->param('save_object_type', 'DnaAlignFeature');
  
  return $runnable;
}

sub results_by_index {
  my ($self, $results) = @_;
  
  $results =~ s/(^Sequence.*\nName.*\n^\-+.*\n)//m;
  my $header = $1;
  
  my @lines = split(/\n/, $results);
  
  my %seqnames;
  foreach my $line (@lines) {
    my ($seqname) = $line =~ /^(\S+)/;
    $seqnames{$seqname}{'result'} .= "$line\n";
  }
  foreach my $seqname (keys %seqnames) {
    $seqnames{$seqname}{'header'} = $header;
  }
  
  return %seqnames;
}

1;
