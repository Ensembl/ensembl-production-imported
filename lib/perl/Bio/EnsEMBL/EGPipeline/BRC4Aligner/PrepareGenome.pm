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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::PrepareGenome;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use File::Spec::Functions qw(catdir);

sub write_output {
  my ($self) = @_;

  my $species = $self->param_required('species');
  my $organism = $self->param_required('organism');
  my $component = $self->param_required('component');
  my $index_dir = $self->param_required('index_dir');
  my $genome_dir = catdir($index_dir, $species);
  my $results_dir = $self->param_required('results_dir');
  my $pipeline_dir = $self->param_required('pipeline_dir');

  my %genome_metadata = (
    species => $species,

    species_results_dir => catdir($results_dir, $component, $organism),

    species_work_dir => catdir($pipeline_dir, $species),
    genome_dir => $genome_dir,
    genome_file => catdir($genome_dir, $species . ".fa"),
    length_file => catdir($genome_dir, $species . ".fa.lengths.txt"),
    genome_bed_file => catdir($genome_dir, $species . ".bed"),
    genome_gff_file => catdir($genome_dir, $species . ".gff"),
    genome_gtf_file => catdir($genome_dir, $species . ".gtf"),
  );

  $self->dataflow_output_id(\%genome_metadata, 1);
}

1;

