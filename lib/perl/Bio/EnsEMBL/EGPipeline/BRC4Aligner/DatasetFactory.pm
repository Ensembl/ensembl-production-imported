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
  my $results_dir = $self->param('results_dir');

  print "Get datasets for $organism from $datasets_file\n";
  my $datasets = $self->get_datasets($datasets_file, $organism);
  print "Number of datasets: " . scalar(@$datasets) . "\n";

  for my $dataset (@$datasets) {
    # Get the a list of run ids from the SRA id
    $dataset->{runs} = $self->all_run_ids($dataset->{runs});

    # Check there is not already a directory for this dataset
    my $ds_dir = catdir($results_dir, $dataset->{name});
    my @files = glob("$ds_dir/*/*");
    
    if (not -e $ds_dir or not @files) {
      $self->dataflow_output_id({ dataset_metadata => $dataset }, 2);
    } else {
      $self->dataflow_output_id({ dataset_metadata => $dataset }, 3);
    }
  }
}

sub all_run_ids {
  my ($self, $runs) = @_;

  for my $run (@$runs) {
    my $accessions = $run->{accessions};

    my @all_run_ids;
    for my $accession (@$accessions) {
      if ($accession =~ /^.RR/) {
        push @all_run_ids, $accession;
      } elsif ($accession =~ /^.RS/) {
        my $adaptor = get_adaptor('Sample');
        for my $sample (@{$adaptor->get_by_accession($accession)}) {
          my @runs = map { $_->accession() } @{$sample->runs()};
          if (not @runs) {
            die("No runs extracted from sample '$accession'");
          }
          push @all_run_ids, @runs;
        }
      } elsif ($accession =~ /^.RP/) {
        my $adaptor = get_adaptor('Study');
        for my $study (@{$adaptor->get_by_accession($accession)}) {
          my @runs = map { $_->accession() } @{$study->runs()};
          if (not @runs) {
            die("No runs extracted from study '$accession'");
          }
          push @all_run_ids, @runs;
        }
      } else {
        die("Unknown SRA accession format: $accession");
      }
    }
    $run->{accessions} = \@all_run_ids;
  }

  return $runs;
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
    } elsif ($dataset->{production_name} eq $organism) {
      push @datasets, $dataset;
    }
  }

  return \@datasets;
}

1;
