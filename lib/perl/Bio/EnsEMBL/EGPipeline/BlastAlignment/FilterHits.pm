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

Bio::EnsEMBL::EGPipeline::BlastAlignment::FilterHits

=head1 DESCRIPTION

Retain the top X hits across a genome, by selecting those with the
best E-value. Note that a single hit might map to multiple locations,
but each of those alignments corresponds to unique protein sequence.
In other words, the 'h' coordinates will not overlap.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::BlastAlignment::FilterHits;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub param_defaults {
  return {
    'db_type'       => 'core',
    'query_type'    => 'dna',
    'database_type' => 'pep',
    'filter_top_x'  => 1,
    'create_gff'    => 0,
  };
}

sub run {
  my ($self) = @_;
  my $db_type       = $self->param_required('db_type');
  my $query_type    = $self->param_required('query_type');
  my $database_type = $self->param_required('database_type');
  my $filter_top_x  = $self->param_required('filter_top_x');
  my $logic_name    = $self->param_required('logic_name');
  
  $self->dbc && $self->dbc->disconnect_if_idle();
  
  if ($filter_top_x) {
    my $dba = $self->get_DBAdaptor($db_type);
    
    if ($database_type eq 'pep') {
      if ($query_type eq 'pep') {
        $self->best_pep_pep($dba, $filter_top_x, $logic_name);
      } else {
        $self->best_pep_dna($dba, $filter_top_x, $logic_name);
      }
    } else {
      $self->best_dna_dna($dba, $filter_top_x, $logic_name);
    }
  }
}

sub write_output {
  my ($self) = @_;
  my $create_gff = $self->param_required('create_gff');
  my $logic_name = $self->param_required('logic_name');
  
  if ($create_gff) {
    $self->dataflow_output_id({'logic_name' => [$logic_name]}, 2);
  }
}

sub best_dna_dna {
  my ($self, $dba, $top_x, $logic_name) = @_;
  
  my %grouped_features;
  
  my $sa = $dba->get_adaptor("Slice");
  my $dafa = $dba->get_adaptor("DNAAlignFeature");
  
  my $slices = $sa->fetch_all('toplevel');
  foreach my $slice (@$slices) {
    my $dafs = $dafa->fetch_all_by_Slice($slice, $logic_name);
    
    $self->process_features($dafs, \%grouped_features);
  }
  
  $self->remove_features($top_x, $dafa, \%grouped_features);
}

sub best_pep_dna {
  my ($self, $dba, $top_x, $logic_name) = @_;
  
  my %grouped_features;
  
  my $sa = $dba->get_adaptor("Slice");
  my $pafa = $dba->get_adaptor("ProteinAlignFeature");
  
  my $slices = $sa->fetch_all('toplevel');
  foreach my $slice (@$slices) {
    my $pafs = $pafa->fetch_all_by_Slice($slice, $logic_name);
    
    $self->process_features($pafs, \%grouped_features);
  }
  
  $self->remove_features($top_x, $pafa, \%grouped_features);
}

sub best_pep_pep {
  my ($self, $dba, $top_x, $logic_name) = @_;
  
  my %grouped_features;
  
  my $pfa = $dba->get_adaptor("ProteinFeature");
  my $pfs = $pfa->fetch_all_by_logic_name($logic_name);
  
  $self->process_features($pfs, \%grouped_features);
  
  $self->remove_features($top_x, $pfa, \%grouped_features);
}

sub process_features {
  my ($self, $features, $grouped_features) = @_;
  
  my %features;
  foreach my $feature (@$features) {
    push @{$features{$feature->hseqname}}, $feature;
  }
  
  # For hits from the same source on a single seq_region, we want
  # to group them into sets within which the hits do not overlap.
  # This means that short hits that would otherwise be lost due to
  # E-value filtering are retained, by piggy-backing on the higher
  # scoring hits.
  foreach my $hit_name (keys %features) {
    my @features = sort {$a->p_value <=> $b->p_value or $b->score <=> $a->score} @{$features{$hit_name}};    
    my $groups = $self->group_features(\@features);
    push @{$$grouped_features{$hit_name}}, @$groups;
  }
}

sub group_features {
  my ($self, $features) = @_;
  my @groups;
  
  while (scalar(@$features)) {
    my $group = [];
    $self->disjoint($features, $features, $group);
    
    push @groups, $group;
  }
  
  return \@groups;
} 

sub disjoint {
  my ($self, $features, $features_subset, $group) = @_;
  
  my $best_hit = shift @$features_subset;
  
  if (scalar(@$group)) {
    my $offset;
    for (my $i=0; $i < scalar(@$features); $i++) {
      if ($$features[$i]->dbID == $best_hit->dbID) {
        $offset = $i;
        last;
      }
    }
    splice(@$features, $offset, 1);
  }
  
  push @$group, $best_hit;
  
  my @disjoint = ();
  if (scalar(@$features_subset)) {
    foreach my $feature (@$features_subset) {
      if (! $self->hit_overlap($best_hit, $feature)) {
        push @disjoint, $feature;
      }
    }
    
    if (scalar(@disjoint)) {
      $self->disjoint($features, \@disjoint, $group);
    }
  }
}

sub remove_features {
  my ($self, $top_x, $adaptor, $grouped_features) = @_;
  
  # Because the features within each group are sorted such that the best
  # hit is first, we can sort the groups by sorting on the first element
  # of each. The top X are then removed from the front of the list of
  # grouped features, and whatever is left hasn't made the cut and is deleted.
  
  foreach my $hit_name (keys %$grouped_features) {
    my @groups = sort {$a->[0]->p_value <=> $b->[0]->p_value or $b->[0]->score <=> $a->[0]->score} @{$$grouped_features{$hit_name}};
    splice(@groups, 0, $top_x);
    foreach my $group (@groups) {
      foreach my $feature (@$group) {
        if ($adaptor->isa('Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor')) {
          my $sth = $adaptor->prepare("DELETE FROM protein_feature WHERE protein_feature_id = ?");
          $sth->execute($feature->dbID);
        } else {
          $adaptor->remove($feature);
        }
      }
    }
  }
}

sub hit_overlap {
  my ($self, $feature1, $feature2) = @_;
  
  # The important thing to note here is that we're checking
  # for overlapping _hit_ coordinates, not seq_region coordinates.
  
  my ($s1, $e1) = ($feature1->hstart, $feature1->hend);
  my ($s2, $e2) = ($feature2->hstart, $feature2->hend);
  
	my $overlap = 0;
  
	if ($s1 >= $s2 && $e1 <= $e2) {
		$overlap = 1;
	} elsif ($s1 <= $s2 && $e1 >= $e2) {
		$overlap = 2;
	} elsif ($e1 >= $s2 && $e1 <= $e2) {
		$overlap = 3;
	} elsif ($s1 >= $s2 && $s1 <= $e2) {
		$overlap = 4;
	}
  
	return $overlap;
}

1;
