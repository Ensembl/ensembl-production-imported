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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::RunFactory;

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
    'tax_id_restrict' => 0,
    'check_library' => 1,
    'default_direction' => 'reverse',
  };
}

sub fetch_input {
  my ($self) = @_;
  
  my $dba = $self->core_dba;
  my $assembly = $dba->get_MetaContainer()->single_value_by_key('assembly.default');
  $self->param('assembly', $assembly);
}

sub run {
  my ($self) = @_;
  
  my $datasets_file = $self->param('datasets_file');
  my $organism = $self->param('organism');
  my $results_dir = $self->param('results_dir');
  my $default_direction = $self->param('default_direction');

  print "Get datasets for $organism from $datasets_file\n";
  my $datasets = $self->get_datasets($datasets_file, $organism);
  print "Number of datasets: " . scalar(@$datasets) . "\n";

  for my $dataset (@$datasets) {
    print "Number of runs for $dataset->{name}: " . scalar(@{ $dataset->{runs} }) . "\n";

    for my $sample (@{ $dataset->{runs} }) {
      
      # Sample data
      my $sample_name = $sample->{name} // $sample->{accessions}->[0];
      die "Missing sample name for $dataset->{name}" if not $sample_name;
      my $sample_dir = catdir($results_dir, $dataset->{name}, $sample_name);
      make_path($sample_dir);

      my $sample_bam_file = catdir($sample_dir, "results.bam");
      (my $sorted_bam_file = $sample_bam_file) =~ s/\.bam/_name.bam/;
      (my $sample_sam_file = $sample_bam_file) =~ s/\.bam/.sam/;

      my $metadata_file = catdir($sample_dir, "metadata.json");

      my @run_ids = $self->runs_from_sra_ids($sample->{accessions});

      my %sample_data = (
        component => $dataset->{component},
        organism => $dataset->{species},
        study_name => $dataset->{name},
        input_is_paired => $sample->{hasPairedEnds} ? 1 : 0,
        input_is_stranded => $sample->{isStrandSpecific} ? 1 : 0,
        trim_reads => $sample->{trim_reads} ? 1 : 0,
        sample_name => $sample_name,
        sample_dir => $sample_dir,
        sample_bam_file => $sample_bam_file,
        sample_sam_file => $sample_sam_file,
        sorted_bam_file => $sorted_bam_file,
        metadata_file => $metadata_file,
      );
      if ($sample_data{input_is_stranded}) {
        if ($sample->{strandDirection}) {
          $sample_data{input_strand_direction} = $sample->{strandDirection};
        } else {
          $sample_data{input_strand_direction} = $default_direction;
        }
        
        if ($sample_data{input_is_paired}) {
          if ($sample_data{input_strand_direction} eq 'forward') {
            $sample_data{input_strandness} = "FR";
          } elsif ($sample_data{input_strand_direction} eq 'reverse') {
            $sample_data{input_strandness} = "RF";
          }
        } else {
          if ($sample_data{input_strand_direction} eq 'forward') {
            $sample_data{input_strandness} = "F";
          } elsif ($sample_data{input_strand_direction} eq 'reverse') {
            $sample_data{input_strandness} = "R";
          }
        }
      }
      
      # In all cases, output the runs metadata (useful for htseq-count)
      $self->dataflow_output_id(\%sample_data, 3);
      
      # Don't remake an existing file
      if (not -s $sorted_bam_file) {
        for my $run_id (@run_ids) {
          my %run_data = (
            run_id => $run_id,
          );
          $self->dataflow_output_id(\%run_data, 2);
        }
        $self->dataflow_output_id(\%sample_data, 4);
      }
    }
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

sub runs_from_sra_ids {
  my ($self, $sra_ids) = @_;
  
  my $tax_id_restrict = $self->param_required('tax_id_restrict');
  
  my @runs;

  print "Number of sra_ids: " . scalar(@$sra_ids) . "\n";
  
  my $run_adaptor = get_adaptor('Run');
  my $study_adaptor = get_adaptor('Study');
  my $exp_adaptor = get_adaptor('Experiment');
  my $sample_adaptor = get_adaptor('Sample');
  foreach my $sra_id (@$sra_ids) {

    if ($sra_id =~ /^[ESD]RP/) {
      # First look for a study
      foreach my $study (@{$study_adaptor->get_by_accession($sra_id)}) {
        foreach my $run (@{$study->runs()}) {
          push @runs, $run;
        }
      }
    }

    if ($sra_id =~ /^[ESD]RX/) {
      # First look for an experiment
      foreach my $exp (@{$exp_adaptor->get_by_accession($sra_id)}) {
        foreach my $run (@{$exp->runs()}) {
          push @runs, $run;
        }
      }
    }

    if ($sra_id =~ /^[ESD]RS/) {
      # First look for a sample
      foreach my $sample (@{$sample_adaptor->get_by_accession($sra_id)}) {
        foreach my $run (@{$sample->runs()}) {
          push @runs, $run;
        }
      }
    }

    # Next look for runs
    if ($sra_id =~ /^[ESD]RR/) {
      foreach my $run (@{$run_adaptor->get_by_accession($sra_id)}) {
        push @runs, $run;
      }
    }

  }
  print "Number of runs: " . scalar(@runs) . "\n";
  
  # Filter to runs to only include the right species and transcriptomic data
  my @filtered_runs;
  foreach my $run (@runs) {
    next if $tax_id_restrict and not $self->tax_id_match($run);
    if ($self->param('check_library') and not $self->is_transcriptomic($run)) {
      print $run->accession . " from @{$sra_ids} is not transcriptomic, skip\n";
      next;
    }
    push @filtered_runs, $run->accession;
  }
  print "Number of filtered runs: " . scalar(@filtered_runs) . "\n";
  
  return @filtered_runs;
}

sub tax_id_match {
  my ($self, $run) = @_;
  
  my $run_tax_id  = $run->sample()->taxon()->taxon_id();
  my $mc          = $self->core_dba->get_MetaContainer();
  my $taxonomy_id = $mc->single_value_by_key('species.taxonomy_id');
  
  return $run_tax_id eq $taxonomy_id;
}

sub is_transcriptomic {
  my ($self, $run) = @_;
  
  my $is_transcriptomic = 0;
  
  # Check study type
  my $study_type = $run->study()->type();
  if ($study_type eq 'Transcriptome Analysis') {
    $is_transcriptomic = 1;
  }
  
  # Otherwise, check experiment type (in case the study is mixed)
  my $design = $run->experiment()->design();
  my $source = $design->{LIBRARY_DESCRIPTOR}->{LIBRARY_SOURCE};
  if ($source eq 'TRANSCRIPTOMIC') {
    $is_transcriptomic = 1;
  }
  
  # Not RNAseq then
  return $is_transcriptomic;
}

1;
