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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::DatasetFactory;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::EnsEMBL::ENA::SRA::BaseSraAdaptor qw(get_adaptor);

use File::Path qw(make_path);
use File::Spec::Functions qw(catdir);
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
  
  my $datasets_file = $self->param('datasets_file');
  my $organism = $self->param('organism');

  print "Get datasets for $organism from $datasets_file\n";
  my $datasets = $self->get_datasets($datasets_file, $organism);
  print "Number of datasets: " . scalar(@$datasets) . "\n";

  for my $dataset (@$datasets) {
    $self->dataflow_output_id({ dataset_metadata => $dataset }, 2);
  }
}

sub get_datasets {
  my ($self, $datasets_file, $organism) = @_;
  
  my $json;
  {
    local $/; #enable slurp
    open my $fh, "<", $datasets_file;
    $json = <$fh>;
  }

  my $data = decode_json($json);

  my @datasets;

  for my $dataset (@$data) {
    if ($dataset->{species} eq $organism) {
      push @datasets, $dataset;
    }
  }

  return \@datasets;
}

1;
