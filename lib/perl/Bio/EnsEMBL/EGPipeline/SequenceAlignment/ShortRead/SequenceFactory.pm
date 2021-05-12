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

package Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::SequenceFactory;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::EnsEMBL::ENA::SRA::BaseSraAdaptor qw(get_adaptor);

use Path::Tiny qw(path);

sub param_defaults {
  my ($self) = @_;
  
  return {
    'reformat_header' => 1,
    'trim_est'        => 1,
    'trimest_exe'     => 'trimest',
    'tax_id_restrict' => 1,
  };
}

sub write_output {
  my ($self) = @_;
  
  my $data_type   = $self->param('data_type');
  my $files       = $self->param('seq_file');
  my $file_pairs  = $self->param('seq_file_pair');
  my $run_ids     = $self->param('run');
  my $study_ids   = $self->param('study');
  my $merge_level = $self->param('merge_level');
  my $merge_id    = $self->param('merge_id');
  my $merge_group = $self->param('merge_group');
  
  $self->param('merge_ids', {});
  
  # From single files
  my @file_output = $self->files($files, $merge_level, $merge_id);
  foreach my $file_output (@file_output) {
    if ($data_type eq 'est') {
      $self->dataflow_output_id($file_output, 5);
    } else {
      $self->dataflow_output_id($file_output, 2);
    }
  }
  
  # From pairs of files
  my @paired_file_output = $self->file_pairs($file_pairs, $merge_level, $merge_id);
  foreach my $paired_file_output (@paired_file_output) {
    $self->dataflow_output_id($paired_file_output, 3);
  }
  
  # From SRA run accessions
  my @run_output = $self->sra_files($run_ids, $study_ids, $merge_level);
  foreach my $run_output (@run_output) {
    $self->dataflow_output_id($run_output, 4);
  }
  
  # From merge groups
  my @merge_group_output = $self->merged_sra_files($merge_group);
  foreach my $merge_group_output (@merge_group_output) {
    $self->dataflow_output_id($merge_group_output, 4);
  }
  
  # Prepare list of tracks to be merged
  my %merge_ids = %{$self->param('merge_ids')};
  $self->throw("No run to align. Abort.") if keys %merge_ids == 0;
  my @output;
  foreach my $merge_id (keys %merge_ids) {
    push @output,
      {
        'merge_id'    => $merge_id,
        'merge_label' => $merge_ids{$merge_id}{'label'},
        'run_ids'     => join(",", @{$merge_ids{$merge_id}{'run_ids'}}),
      };
  }
  
  $self->dataflow_output_id({'merges' => \@output}, 1);
}

sub files {
  my ($self, $files, $merge_level, $merge_id) = @_;
  
  my $data_type       = $self->param('data_type');
  my $reformat_header = $self->param('reformat_header');
  my $trim_est        = $self->param('trim_est');
  my $trimest_exe     = $self->param('trimest_exe');
  
  my @all;
  foreach my $file (@$files) {
    push @all, path($file)->basename;
  }
  
  my @output = ();
  
  foreach my $file (@$files) {
    my $file_varname = 'seq_file_1';
    
    if ($data_type eq 'est') {
      $file_varname = 'fasta_file';
      
      if ($trim_est) {
        $file = $self->trim_est($file);
      }
      if ($reformat_header) {
        $self->reformat_header($file);
      }
    }
    
    if (lc($merge_level) eq 'taxon') {
      if (!defined $merge_id) {
        $merge_id = join('_', @all);
      }
    } else {
      $merge_id = path($file)->basename;
    }
    
    $self->param('merge_ids')->{$merge_id}{'label'} = $merge_id;
    $self->param('merge_ids')->{$merge_id}{'run_ids'} = [];
    
    push @output,
      {
        $file_varname => $file,
        'merge_id'    => $merge_id,
      };
  }
  
  return @output;
}

sub file_pairs {
  my ($self, $file_pairs, $merge_level, $merge_id) = @_;
  
  my @all;
  foreach my $file_pair (@$file_pairs) {
    my @file_pair = split(/\s*,\s*/, $file_pair);
    if (scalar(@file_pair) == 2) {
      (my $file_name = path($file_pair[0])->basename) =~ s/_1\.\w+$//;
      push @all, $file_name;
    } else {
      $self->throw("File pair does not contain a comma-separated pair of files: '$file_pair'");
    }
  }
  
  my @output = ();
  
  foreach my $file_pair (@$file_pairs) {
    my ($seq_file_1, $seq_file_2) = split(/\s*,\s*/, $file_pair);
    
    if (lc($merge_level) eq 'taxon') {
      if (!defined $merge_id) {
        $merge_id = join('_', @all);
      }
    } else {
      ($merge_id = path($seq_file_1)->basename) =~ s/_1\.\w+$//;
    }
    
    $self->param('merge_ids')->{$merge_id}{'label'} = $merge_id;
    $self->param('merge_ids')->{$merge_id}{'run_ids'} = [];
    
    push @output,
      {
        'seq_file_1' => $seq_file_1,
        'seq_file_2' => $seq_file_2,
        'merge_id'   => $merge_id,
      };
  }
  
  return @output;
}

sub sra_files {
  my ($self, $run_ids, $study_ids, $merge_level) = @_;
  
  my $tax_id_restrict = $self->param_required('tax_id_restrict');
  
  my @runs;
  
  my $run_adaptor = get_adaptor('Run');
  foreach my $run_id (@$run_ids) {
    foreach my $run (@{$run_adaptor->get_by_accession($run_id)}) {
      push @runs, $run;
    }
  }
  
  my $study_adaptor = get_adaptor('Study');
  foreach my $study_id (@$study_ids) {
    foreach my $study (@{$study_adaptor->get_by_accession($study_id)}) {
      foreach my $run (@{$study->runs()}) {
        push @runs, $run;
      }
    }
  }
  
  my @output = ();
  
  foreach my $run (@runs) {
    if ($tax_id_restrict) {
      next if not $self->tax_id_match($run);
    }
    next if not $self->is_transcriptomic($run);
    
    my $run_id = $run->accession;
    my ($merge_id, $merge_label) = $self->sra_merge_id($merge_level, $run);
    
    $self->param('merge_ids')->{$merge_id}{'label'} = $merge_label;
    push @{$self->param('merge_ids')->{$merge_id}{'run_ids'}}, $run_id;
    
    push @output,
      {
        'run_id'   => $run_id,
        'merge_id' => $merge_id,
      };
  }
  
  return @output;
}

sub merged_sra_files {
  my ($self, $merge_groups) = @_;
  
  my $tax_id_restrict = $self->param_required('tax_id_restrict');
  my $run_adaptor = get_adaptor('Run');
  my @output = ();
  
  foreach my $merge_id (keys %$merge_groups) {
      my @run_ids = split/,/, $merge_groups->{$merge_id};
      my @runs;
      
      # Get run object from SRA accession
      foreach my $run_id (@run_ids) {
          foreach my $run (@{$run_adaptor->get_by_accession($run_id)}) {
              push @runs, $run;
          }
      }
      
      # Add each run to the group
      foreach my $run (@runs) {
          if ($tax_id_restrict) {
              next if not $self->tax_id_match($run);
          }
          next if not $self->is_transcriptomic($run);
          
          my $run_id = $run->accession;
          my $merge_label = $merge_id;

          $self->param('merge_ids')->{$merge_id}{'label'} = $merge_label;
          
          # Store the list of runs for this merge_id (for the MergeFactory)
          push @{$self->param('merge_ids')->{$merge_id}{'run_ids'}}, $run_id;
            
          # And add this run to be processed
          push @output,
          {
              'run_id'   => $run_id,
              'merge_id' => $merge_id,
          };
      }
  }
  
  return @output;
}

sub reformat_header {
  my ($self, $fasta_file) = @_;
  
  my $file = path($fasta_file);
  my $data = $file->slurp;
  $data =~ s/^>gi\|\d+\|gb\|([^\|]+)\|\S+/>$1/g;
  $file->spew($data);
}

sub trim_est {
  my ($self, $fasta_file) = @_;
  
  my $trimest_exe = $self->param('trimest_exe');
  (my $trimmed_file = $fasta_file) =~ s/(\.\w+)$/-trimmed$1/;
  
  my $cmd = "$trimest_exe -seq $fasta_file -out $trimmed_file";
  system($cmd) == 0 || $self->throw("Cannot execute $cmd");
  
  return $trimmed_file;
}

sub sra_merge_id {
  my ($self, $merge_level, $run) = @_;
  
  my ($merge_id, $merge_label);
  
  if ($merge_level =~ /study/i) {
    $merge_id = $run->study->accession;
    $merge_label = $run->study->title;
    
  } elsif ($merge_level =~ /taxon/i) {
    my $taxon = $run->sample->taxon;
    $merge_id = $self->param('merge_id');
    
    if ($taxon) {
      $merge_id //= $taxon->taxon_id;
      my $taxonomy_adaptor = Bio::EnsEMBL::ENA::SRA::BaseSraAdaptor->new->taxonomy_adaptor();
      my $node = $taxonomy_adaptor->fetch_by_taxon_id($taxon->taxon_id);
      $merge_label = $node->name;
    } else {
      $merge_id //= 'unknown';
      $merge_label = 'unknown species';
    }
    
  } elsif ($merge_level =~ /sample/i) {
    $merge_id = $run->sample->accession;
    $merge_label = $run->sample->title // $run->sample->description // '';
    
  } elsif ($merge_level =~ /experiment/i) {
    $merge_id = $run->experiment->accession;
    $merge_label = $run->experiment->title;
    
  } elsif ($merge_level =~ /run/i) {
    $merge_id = $run->accession;
    $merge_label = $run->title;
    
  } else {
    $self->throw("Unrecognised merge level, '$merge_level'");
  }
  
  return ($merge_id, $merge_label);
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
