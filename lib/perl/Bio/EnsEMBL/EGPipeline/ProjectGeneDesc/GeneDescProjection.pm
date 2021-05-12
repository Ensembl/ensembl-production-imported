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


=pod

=head1 NAME

Bio::EnsEMBL::EGPipeline::ProjectGeneDesc::GeneDescProjection

=head1 DESCRIPTION

Project gene descriptions from one species to another using the orthologies 
from the Compara ProteinTree pipeline. 

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::ProjectGeneDesc::GeneDescProjection;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub run {
  my ($self) = @_;
  
  my $from_species = $self->param_required('source');
  my $to_species   = $self->param_required('species');
  my $log_file     = $self->param_required('log_file');
  my $store_data   = $self->param_required('store_data');
  
  my $from_ga = Bio::EnsEMBL::Registry->get_adaptor($from_species, 'core', 'Gene');
  my $to_ga   = Bio::EnsEMBL::Registry->get_adaptor($to_species, 'core', 'Gene');
  $self->throw("Cannot connect to $from_species core database") if !$from_ga;
  $self->throw("Cannot connect to $to_species core database") if !$to_ga ;
  
  my $from_db    = $from_ga->dbc->dbname;
  my $to_db      = $to_ga->dbc->dbname;
  my $from_count = $from_ga->count_all_by_biotype('protein_coding');
  my $to_count   = $to_ga->count_all_by_biotype('protein_coding');
  
  open my $data, '>', $log_file or $self->throw($!);
  print $data "\nProjection log\n";
  print $data "\tfrom: $from_species ($from_db), $from_count protein-coding genes\n";
  print $data "\tto:   $to_species ($to_db), $to_count protein-coding genes\n";
  
  if ($store_data) {
    $self->remove_existing($from_species);
  }
  
  my $homologies = $self->fetch_homologies($from_species, $to_species);
  my $homology_count = scalar keys(%$homologies);
  print $data "Fetched $homology_count homologies\n";
  
  my $projected = 0;
  
  foreach my $from_stable_id (keys %$homologies) {
    my $from_gene = $from_ga->fetch_by_stable_id($from_stable_id);
    if ($from_gene) {
      my @to_genes = @{$$homologies{$from_stable_id}};
      
      foreach my $to_stable_id (@to_genes) {
        my $to_gene = $to_ga->fetch_by_stable_id($to_stable_id);
        if ($to_gene) {
          $projected += $self->project_description($from_gene, $to_gene, $to_ga, $data);
        }
      }
    }
  }
  
  close $data;
  
  my $summary =
    "From $from_species ($from_count genes) ".
    "to $to_species ($to_count genes): ".
    "$homology_count homologies; ".
    "$projected projected descriptions";
  $self->warning($summary);
  
  my @summary = (
    $from_species,
    $from_count,
    $to_species,
    $to_count,
    $homology_count,
    $projected,
  );
  
  $self->param('summary', \@summary)
}

sub write_output {
  my ($self)  = @_;
  
  $self->dataflow_output_id({'summary' => $self->param('summary')}, 1);
}

sub fetch_homologies {
  my ($self, $from_species, $to_species) = @_;
  
  my $compara_name    = $self->param('compara_name');
  my $compara_db_name = $self->param('compara_db_name');
  my $method_link     = $self->param_required('method_link');
  
  my ($gdba, $mlssa, $ha);
  if ($compara_db_name) {
    my $dbas = Bio::EnsEMBL::Registry->get_all_DBAdaptors_by_dbname($compara_db_name);
    foreach my $dba (@$dbas) {
      $gdba  = $dba->get_adaptor('GenomeDB');
      $mlssa = $dba->get_adaptor('MethodLinkSpeciesSet');
      $ha    = $dba->get_adaptor('Homology');
      last;
    }
  } else {
    if (!$compara_name) {
      $self->throw("Either -compara_name or -compara_db_name must be given");
    }
    $gdba  = Bio::EnsEMBL::Registry->get_adaptor($compara_name, 'compara', 'GenomeDB');
    $mlssa = Bio::EnsEMBL::Registry->get_adaptor($compara_name, 'compara', 'MethodLinkSpeciesSet');
    $ha    = Bio::EnsEMBL::Registry->get_adaptor($compara_name, 'compara', 'Homology');
  }
  $self->throw("Cannot connect to $compara_name compara database") if !$gdba;
  
  my $from_gdb   = $gdba->fetch_by_registry_name($from_species);
  my $to_gdb     = $gdba->fetch_by_registry_name($to_species);
  my $mlss       = $mlssa->fetch_by_method_link_type_GenomeDBs($method_link, [$from_gdb, $to_gdb]);
  my $homologies = $ha->fetch_all_by_MethodLinkSpeciesSet($mlss);
  
  return $self->filter_homologies($from_species, $homologies);
}

sub filter_homologies {
  my ($self, $from_species, $homologies) = @_;

  my $homology_types     = $self->param_required('homology_types');
  my $percent_id_filter  = $self->param_required('percent_id_filter');
  my $percent_cov_filter = $self->param_required('percent_cov_filter');
  
  my %homology_types = map {$_ => 1} @$homology_types;
  
  my %filtered_homologies;
  foreach my $homology (@{$homologies}) {
    if (exists $homology_types{$homology->description}) {
      next unless $homology->is_tree_compliant() == 1;
      
      my ($from_stable_id, @to_stable_ids);
      my $from_seen = 0;
      my $members = $homology->get_all_Members();
      foreach my $member (@$members) {
        if ($member->perc_id() >= $percent_id_filter) {
          if ($member->perc_cov() >= $percent_cov_filter) {
            my $gene_member = $member->gene_member();
            if ($gene_member->genome_db->name() eq $from_species) {
              $from_stable_id = $gene_member->stable_id();
              $from_seen++;
            } else {
              push @to_stable_ids, $gene_member->stable_id();
            }
          }
          if ($from_seen == 1) {
            $filtered_homologies{$from_stable_id} = \@to_stable_ids;
          }
        }
      }
    }
  }
  
  return \%filtered_homologies;
}

sub remove_existing {
  my ($self, $from_species) = @_;
  
  my $from_species_text = ucfirst($from_species);
  $from_species_text =~ s/_/ /g;
  
  my $dbh = $self->core_dbh();
  my $sql =
    "UPDATE gene SET description = NULL ".
    "WHERE description LIKE '%[Source:Projected from $from_species_text%'";
  
  $dbh->do($sql) or $self->throw("Failed to execute: $sql");
}

sub project_description {
  my ($self, $from_gene, $to_gene, $to_ga, $data) = @_;
  
  my $store_data     = $self->param_required('store_data');
  my $ignore_source  = $self->param_required('ignore_source');
  my $replace_target = $self->param_required('replace_target');
  
  my $projected = 0;
  
  my $from_desc = $from_gene->description();  
  if ($from_desc) {
    my $ignore = grep {$from_desc =~ /\Q$_\E/} @$ignore_source;
    if (!$ignore) {
      my $to_desc = $to_gene->description();
      my $replace = grep {$to_desc =~ /\Q$_\E/} @$replace_target if $to_desc;
      if (!$to_desc || $replace) {
        my $from_stable_id = $from_gene->stable_id();
        my $to_stable_id = $to_gene->stable_id();
        my $from_species_text = ucfirst($from_gene->species());
        $from_species_text =~ s/_/ /g;
        $from_desc = $self->scrub_description($from_gene, $from_desc);
        $from_desc =~ s/(\[Source:)/$1Projected from $from_species_text ($from_stable_id) /;
        
        $projected = 1;
        print $data "\t$from_stable_id\t$to_stable_id\t$from_desc\n";
        
        if ($store_data) {
          $to_gene->description($from_desc);
          $to_ga->update($to_gene);
        }
      }
    }
  }
  
  return $projected;
}

sub scrub_description {
  my ($self, $from_gene, $from_desc) = @_;
  
  # Remove gene names from descriptions
  my $gene_name = $from_gene->external_name;
  if (defined $gene_name) {
    $from_desc =~ s/[,\s]*\($gene_name\)//;
  }
  
  # Remove isolated digit suffixes
  $from_desc =~ s/\s+\d+(\s+\[Source:)/$1/;
  
  return $from_desc;
}

1
