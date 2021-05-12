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

package Bio::EnsEMBL::EGPipeline::RNAFeatures::TaxonomicFilter;

use strict;
use warnings;

use Path::Tiny qw(path);
use File::Basename;

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub param_defaults {
  return {
    'rfam_rrna'           => 1,
    'rfam_trna'           => 0,
    'rfam_blacklist'      => [],
    'taxonomic_filtering' => 1,
    'taxonomic_lca'       => 0,
    'taxonomic_levels'    => [],
    'taxonomic_threshold' => 0.02,
    'taxonomic_minimum'   => 50,
  };
}

sub run {
  my ($self) = @_;
  
  my $rfam_cm_file       = $self->param_required('rfam_cm_file');
  my $filtered_cm_file   = $self->param_required('filtered_cm_file');
  my $rfam_rrna          = $self->param('rfam_rrna');
  my $rfam_trna          = $self->param('rfam_trna');
  my $rfam_blacklist     = $self->param('rfam_blacklist');
  my $rfam_whitelist     = $self->param('rfam_whitelist');
  my $rfam_taxonomy_file = $self->param('rfam_taxonomy_file');
  my $filtering          = $self->param('taxonomic_filtering');
  my $lca                = $self->param('taxonomic_lca');
  my $levels             = $self->param('taxonomic_levels');
  my $threshold          = $self->param('taxonomic_threshold');
  my $minimum            = $self->param('taxonomic_minimum');
  
  my $dba = $self->core_dba();
  
  my $rfam_species_path = dirname $filtered_cm_file;
  path($rfam_species_path)->mkpath;
  path($rfam_cm_file)->copy($filtered_cm_file);
  
  # Remove any existing indexes
  unlink glob "$filtered_cm_file.i*";
  
  my $cm_path = path($filtered_cm_file);
  my $cm = $cm_path->slurp;
  
  if (!$rfam_rrna) {
    $cm = $self->filter_rrna($cm);
  }
  
  if (!$rfam_trna) {
    $cm = $self->filter_trna($cm);
  }
  
  $cm = $self->filter_blacklist($cm, $rfam_blacklist, $dba);
  
  if ($filtering) {
    $cm = $self->filter_levels(
      $cm, $rfam_taxonomy_file, $rfam_whitelist, $lca, $levels, $threshold, $minimum, $dba);
  }
  
  my $db_connection = $dba->dbc();
  $db_connection && $db_connection->disconnect_if_idle();

  $cm_path->spew($cm);
}

sub filter_rfam {
  my ($self, $cm, $rfam_acc) = @_;
  
  my @filtered;
  my @cm_models = split(/\/\/\n/, $cm);
  
  my %rfam_acc = map { $_ => 1 } @$rfam_acc;
  foreach my $cm_model (@cm_models) {
    my ($rfam_acc) = $cm_model =~ /^ACC\s+(\S+)/m;
    
    if (! exists($rfam_acc{$rfam_acc})) {
      push @filtered, $cm_model;
    } else {
      $self->warning("Rfam model $rfam_acc removed by filtering.");
    }
  }
  
  return join("//\n", @filtered)."//\n";
}

sub filter_rrna {
  my ($self, $cm) = @_;
  my @rfam_acc = $cm =~ /^ACC\s+(\S+)\nDESC\s+.*ribosomal RNA.*/gm;
  return $self->filter_rfam($cm, \@rfam_acc);
}

sub filter_trna {
  my ($self, $cm) = @_;
  my @rfam_acc = $cm =~ /^NAME\s+tRNA.*\nACC\s+(\S+)/gm;
  return $self->filter_rfam($cm, \@rfam_acc);
}

sub filter_blacklist {
  my ($self, $cm, $blacklist, $dba) = @_;
  my @rfam_acc;
  foreach my $name (keys %$blacklist) {
    if ($self->has_ancestor($dba, $name)) {
      push @rfam_acc, @{$$blacklist{$name}};
    }
  }
  return $self->filter_rfam($cm, \@rfam_acc);
}

sub filter_levels {
  my ($self, $cm, $file, $whitelist, $lca, $levels, $threshold, $minimum, $dba) = @_;
  
  my ($rows, $columns) = $self->parse_rfam_taxonomy_file($file);
  
  if (scalar(@$levels) == 0) {
    my $division = $dba->get_MetaContainer->get_division() || 'Ensembl';
    push @$levels, $division;
  }
  
  my %whitelist;
  foreach my $level (keys %$whitelist) {
    foreach my $rfam_acc (@{$$whitelist{$level}}) {
      $whitelist{$rfam_acc}{$level} = 1;
    }
  }
  
  my @rfam_acc;
  foreach my $row (@$rows) {
    my @row = split(/\t/, $row);
    my $rfam_acc = $row[$$columns{'Rfam_acc'}];
    my $total;
    if (exists $$columns{'Taxa'}) {
      $total = $row[$$columns{'Taxa'}];
    } else {
      $total = $row[$$columns{'Sequences'}];
    }
    my $valid = 0;
    
    if ($total > 0) {
      foreach my $level (@$levels) {
        my $level_total = $row[$$columns{$level}];
        
        if (exists($whitelist{$rfam_acc}{$level})) {
          $valid = 1;
          last;
        } elsif (($level_total/$total >= $threshold) || ($level_total > $minimum)) {
          if ($lca) {
            my $lca_level = $row[$$columns{"LCA_$level"}];
            if ($self->has_ancestor($dba, $lca_level)) {
              $valid = 1;
              last;
            }
          } else {
            $valid = 1;
            last;
          }
        }
      }
    }
    
    if (!$valid) {
      push @rfam_acc, $rfam_acc;
    }
  }
  
  return $self->filter_rfam($cm, \@rfam_acc);
}

sub parse_rfam_taxonomy_file {
  my ($self, $rfam_taxonomy_file) = @_;
  
  my $rfam_taxonomy_path = path($rfam_taxonomy_file);
  my $rfam_taxonomy = $rfam_taxonomy_path->slurp;
  
  my @rows = split(/\n/, $rfam_taxonomy);
  my $header = shift @rows;
  my @columns = split(/\t/, $header);
  my $cursor = 0;
  my %columns = map { $_ => $cursor++ } @columns;
  
  return (\@rows, \%columns);
}

sub has_ancestor {
  my ($self, $dba, $name) = @_;
  my $meta_key = 'species.classification';
  my $meta_values = $dba->get_MetaContainer->list_value_by_key($meta_key);
  
  if (exists {map {$_ => 1} @$meta_values}->{$name}) {
    return 1;
  }
  return 0;
}

1;
