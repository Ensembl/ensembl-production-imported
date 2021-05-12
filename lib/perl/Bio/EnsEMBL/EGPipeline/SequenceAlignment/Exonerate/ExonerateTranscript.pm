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

package Bio::EnsEMBL::EGPipeline::SequenceAlignment::Exonerate::ExonerateTranscript;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript;
use Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter;

use base ('Bio::EnsEMBL::EGPipeline::SequenceAlignment::Exonerate::Exonerate');

sub param_defaults {
  my $self = shift @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'coverage'       => 90,
    'percent_id'     => 97,
    'best_in_genome' => 1,
  };
}

sub fetch_runnable {
  my $self = shift @_;
  
  my %parameters;
  if (%{$self->param('parameters_hash')}) {
    %parameters = %{$self->param('parameters_hash')};
  }
  
  my $use_exonerate_server = $self->param('use_exonerate_server'),
  my $target_file = $self->param_required('queryfile');
  my $server_file = $self->param_required('server_file');
  
  if ($use_exonerate_server) {
    open my $fh, $server_file or $self->throw("Failed to open $server_file");
    my $server = <$fh>;
    close $fh;
    $parameters{'TARGET_FILE'} = $server;
  } else {
    $parameters{'TARGET_FILE'} = $target_file;
  }
  
  $parameters{'QUERY_FILE'} = $self->param_required('seq_file');
  $parameters{'QUERY_TYPE'} = $self->param_required('seq_type');
  $parameters{'COVERAGE_BY_ALIGNED'} = 1;
  
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript->new
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
  
  $self->param('parse_filehandle', 1);
  $self->param('output_not_set', 1);
  
  # Add a run_analysis sub, since the module lacks one and the generic
  # one from Runnable.pm doesn't quite work.
  sub Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript::run_analysis {
    my ($self, $program) = @_;
    
    if (!$program) {
      $program = $self->program;
    }
    $self->throw($program." is not executable.") unless($program && -x $program);
    
    $self->resultsfile($self->query_file.'.out');
    my $err_file = $self->query_file.'.err';
    
    my $command = $program." ";
    $command .= $self->options." " if($self->options);
    $command .= " --querytype ".  $self->query_type;
    $command .= " --targettype ". $self->target_type;
    $command .= " --query ".      $self->query_file;
    $command .= " --target ".     $self->target_file;
    $command .= " > ".            $self->resultsfile;
    $command .= " 2> ".           $err_file;
    
    if (!-e $self->resultsfile) {
      system("exec $command") == 0 or $self->throw("Failed to run $command, errors in $err_file");
    }
  }
  
  $self->param('save_object_type', 'Gene');
  
  return $runnable;
}

sub filter_output {
  my ($self, $runnable) = @_;
  
  my $filter = Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter->new
    (
      -coverage                 => $self->param('coverage'),
      -percent_id               => $self->param('percent_id'),
      -best_in_genome           => $self->param('best_in_genome'),
      -reject_processed_pseudos => 1,
    );
  
  my $filtered = $filter->filter_results($runnable->output);
  
  # Output is cumulative, so need to manually erase existing results.
  # (Note that calling the runnable's 'output' method will NOT work.
  $runnable->{'output'} = $filtered;
  
  # After filtering, generate genes for saving to the db.
  $runnable->{'output'} = $self->make_genes($runnable);
}

sub make_genes {
  my ($self, $runnable) = @_;
  
  my $slice = $runnable->query;
  my $analysis = $self->param('analysis');
  my $biotype = $self->param('biotype');
  my @genes;
  
  foreach my $transcript (@{$runnable->output}) {
    $transcript->slice($slice);
    $transcript->analysis($analysis);
    $transcript->biotype($biotype);
    
    foreach my $sf (@{$transcript->get_all_supporting_features}) {
      $sf->slice($slice);
      $sf->analysis($analysis);
    }
    
    foreach my $exon (@{$transcript->get_all_Exons}) {
      $exon->slice($slice);
      foreach my $sf (@{$exon->get_all_supporting_features}) {
        $sf->slice($slice);
        $sf->analysis($analysis);
      }
    }
    
    my $gene = Bio::EnsEMBL::Gene->new(
      -analysis => $analysis,
      -biotype  => $biotype,
    );
    $gene->add_Transcript($transcript);
    push @genes, $gene;
  }
  
  return \@genes;
}

1;
