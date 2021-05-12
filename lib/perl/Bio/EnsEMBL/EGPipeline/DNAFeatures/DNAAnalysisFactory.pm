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

package Bio::EnsEMBL::EGPipeline::DNAFeatures::DNAAnalysisFactory;

use strict;
use warnings;
use Bio::EnsEMBL::EGPipeline::DNAFeatures::Utils qw(get_repbase_species_lib);
use File::Spec::Functions qw(catfile);

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub param_defaults {
  return {
    'dust'              => 1,
    'repeatmasker'      => 1,
    'redatrepeatmasker' => 1,
    'trf'               => 1,
    'max_seq_length'    => 1000000,
    'species_list'      => [],
    'guess_repbase_species' => 0,
  };
}

sub run {
  my ($self) = @_;
  my $species = $self->param('species');
  my $species_list = $self->param_required('species_list');
  my $dna_analyses = $self->param_required('dna_analyses');
  my $always_use_repbase = $self->param_required('always_use_repbase');
  my $rd_rm_library = $self->param_required('rd_rm_library');
  my $rm_library = $self->param_required('rm_library');
  my $rm_sensitivity = $self->param_required('rm_sensitivity');
  my $rm_logic_name = $self->param_required('rm_logic_name');
  my $pipeline_dir = $self->param_required('pipeline_dir');
  my $db_backup_file = $self->param('db_backup_file');
  my $guess_repbase_species = $self->param('guess_repbase_species');
  
  # More than one species implies a collection db, and we cannot 
  # do anything species-specific.
  my $species_specific = 1;
  if (scalar(@$species_list) > 1 || ! defined $species) {
    $species_specific = 0;
  }
  
  my $filtered_analyses = [];

  foreach my $analysis (@{$dna_analyses}) {

    my $logic_name = $$analysis{'logic_name'};

    if ($logic_name eq 'dust') {
        push @$filtered_analyses, $analysis
          if $self->param('dust');
    }

    elsif ($logic_name eq 'trf') {
        push @$filtered_analyses, $analysis
          if $self->param('trf');
    }

    elsif ($logic_name eq 'repeatmask_repbase') {
      if ($self->param('repeatmasker')) {
        # This is analysis with the repbase library. If the
        # sensitivity is given then use that setting, otherwise the
        # default is 'medium'.
        my $sensitivity = 'automatic';
        
        if ($species_specific) {

          my $repbase_species;
          
          if ($guess_repbase_species) {
            my $temp_lib_file = catfile($pipeline_dir, $species . ".repbase.lib.tmp");
            $repbase_species = get_repbase_species_lib($self->core_dba(), $temp_lib_file);
          } else {
            # use $species name if not overriden by repbase_species
            my $repbase_species =
              ( exists $$analysis{'repbase_species'} && $$analysis{'repbase_species'} )
                ? $$analysis{'repbase_species'}
                : $species;

            ##Special case for strains with longer names (e.g. triticum_aestivum_stanley)
            if ($repbase_species =~ /(\w+_\w+)_/){
                $repbase_species =$1;
            }
            
            # RepeatMasker requires a species name with a space.
            $repbase_species = ucfirst($repbase_species);
            $repbase_species =~ s/_/ /;
          }

          if ($always_use_repbase || (! exists $$rm_library{$species} && ! exists $$rm_library{'all'})) {
            if ($self->check_repeatmasker($analysis, $repbase_species, $pipeline_dir)) {
              $$analysis{'parameters'} .= " -species \"$repbase_species\" ";
            }

            if (exists $$rm_sensitivity{$species} || exists $$rm_sensitivity{'all'}) {
              $sensitivity = $$rm_sensitivity{$species} || $$rm_sensitivity{'all'};
            }
          }
        }
        else {
          if ($always_use_repbase || (! exists $$rm_library{'all'})) {
            if (exists $$rm_sensitivity{'all'}) {
              $sensitivity = $$rm_sensitivity{'all'};
            }
          }
        }
        
        if ($sensitivity =~ /^automatic$/i) {
          $sensitivity = 'medium';
        }
        $self->warning("Sensitivity set to $sensitivity for RepBase library");
        $$analysis{'parameters'} .= $self->rm_engine_params($sensitivity);
          
        push @$filtered_analyses, $analysis;
      }
    }

    elsif ($logic_name eq 'repeatmask_redat') {
      if ($self->param('redatrepeatmasker')) {
        $$analysis{'db_file'} = $rd_rm_library;
        $$analysis{'parameters'} .= " -lib \"$rd_rm_library\" ";
        my $sensitivity = 'automatic';
        
        if (exists $$rm_sensitivity{$species} || exists $$rm_sensitivity{'all'}) {
          $sensitivity = $$rm_sensitivity{$species} || $$rm_sensitivity{'all'};
        }
        else {
          $sensitivity = 'low';
        }
        $self->warning("Sensitivity for $species set to $sensitivity for the REdat library");
        $$analysis{'parameters'} .= $self->rm_engine_params($sensitivity);
          
        push @$filtered_analyses, $analysis;
      }
    }

    elsif ($logic_name eq 'repeatmask_customlib') {
      if ($self->param('repeatmasker')) {
        # This is analysis with a custom repeat library. Such a
        # library can be species-specific, or a single library can be
        # used against all species to which the pipeline is applied.
        # Similarly, the logic_name for the analysis and the
        # sensitivity can be species-specific or universal. If the
        # sensitivity is set to 'automatic' then the file sizes are
        # used to have a stab at suitable parameters; note that this
        # is species-specific, even if a single repeat library is
        # applied to multiple species.
        
        if (exists $$rm_library{$species} || exists $$rm_library{'all'}) {
          my $library = $$rm_library{$species} || $$rm_library{'all'};
          $$analysis{'db_file'} = $library;
          $$analysis{'parameters'} .= " -lib \"$library\" ";
          
          my $sensitivity = 'automatic';
          if (exists $$rm_sensitivity{$species} || exists $$rm_sensitivity{'all'}) {
            $sensitivity = $$rm_sensitivity{$species} || $$rm_sensitivity{'all'};
          }
          else {
            $sensitivity = $self->set_sensitivity($library);
          }
          $self->warning("Sensitivity for $species set to $sensitivity for library '$library'");
          $$analysis{'parameters'} .= $self->rm_engine_params($sensitivity);
          
          if (exists $$rm_logic_name{$species} || exists $$rm_logic_name{'all'}) {
            my $name = $$rm_logic_name{$species} || $$rm_logic_name{'all'};
            $$analysis{'logic_name'} = $name;
          }
          
          push @$filtered_analyses, $analysis;
        }
      }
    }

    else {
      $self->warning("Unrecognised analysis '$logic_name' in 'dna_analyses' parameter.");
    }
  }
  
  $self->core_dbc()->disconnect_if_idle();
  
  foreach my $analysis (@{$filtered_analyses}) {
    $$analysis{'species'} = $species;
    $$analysis{'db_backup_file'} = $db_backup_file if defined $db_backup_file;
  }
  
  $self->param('filtered_analyses', $filtered_analyses);
  
}

sub write_output {
  my ($self) = @_;
  
  $self->dataflow_output_id($self->param('filtered_analyses'), 2);
  
}

sub repeatmasker_species {
  
}

sub repeatmasker_generic {
  
}

# If RepeatMasker hasn't got a species in its taxonomy, it will fail.
sub check_repeatmasker {
  my ($self, $analysis, $species_rm, $pipeline_dir) = @_;
  
  my $program_file = $$analysis{'program_file'};
  my $parameters = $$analysis{'parameters'};
  $parameters .= " -species \"$species_rm\"";
  
  my $test_fasta_file = "$pipeline_dir/rm_test.fa";
  open (FASTA, ">$test_fasta_file") or die "Failed to create test file '$test_fasta_file'";
  print FASTA ">test\nATTATT\n";
  close (FASTA); 
  my $rm_out = `$program_file $parameters $test_fasta_file 2>&1`;
  unlink $test_fasta_file;
  
  if ($rm_out =~ /Species.*is not known to RepeatMasker/) {
    return 0;
  } else {
    return 1;
  }
}

sub rm_engine_params {
  my ($self, $sensitivity) = @_;
  my $rm_engine_params = '';
  
  if ($sensitivity =~ /^very_?low$/i) {
    $rm_engine_params = ' -engine ncbi -q ';
  } elsif ($sensitivity =~ /^low$/i) {
    $rm_engine_params = ' -engine ncbi ';
  } elsif ($sensitivity =~ /^medium$/i) {
    $rm_engine_params = ' -engine crossmatch -q ';
  } elsif ($sensitivity =~ /^high$/i) {
    $rm_engine_params = ' -engine crossmatch ';
  } elsif ($sensitivity =~ /^very_?high$/i) {
    $rm_engine_params = ' -engine crossmatch -s ';
  }
  
  return $rm_engine_params;
}

sub set_sensitivity {
  my ($self, $library) = @_;
  my $large_genome = 1e9;
  my $large_genome_medium_limit = 50e12;
  my $small_genome_high_limit = 15e12;
  my $small_genome_medium_limit = 100e12;
  
  my $sensitivity;
  
  # Get genome size.
  my $dbh = $self->core_dbh();
  my $seq_region_length_sql = '
    SELECT sum(length) FROM
      seq_region INNER JOIN
      seq_region_attrib USING (seq_region_id) INNER JOIN
      attrib_type USING (attrib_type_id)
    WHERE code = "toplevel"
  ;';
  my ($seq_region_length) = $dbh->selectrow_array($seq_region_length_sql);
  
  # Get product of input and library size, as a rough guide
  # to how much 'stuff' will need to be done by RepeatMasker.
  my $library_size = -s $library;
  my $size = $self->param('max_seq_length') * $library_size;
  
  if ($seq_region_length > $large_genome) {
    if ($size < $large_genome_medium_limit) {
      $sensitivity = 'medium';
    } else {
      $sensitivity = 'low';
    }
  } else {
    if ($size < $small_genome_high_limit) {
      $sensitivity = 'high';
    } elsif ($size < $small_genome_medium_limit) {
      $sensitivity = 'medium';
    } else {
      $sensitivity = 'low';
    }
  }
  
  return $sensitivity;
}

1;
