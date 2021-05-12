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

Bio::EnsEMBL::EGPipeline::RNAFeatures::CreateCmscanGenes

=cut

package Bio::EnsEMBL::EGPipeline::RNAFeatures::CreateCmscanGenes;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::RNAFeatures::CreateGenes');

sub param_defaults {
  my ($self) = @_;
  return {
    %{$self->SUPER::param_defaults},
    'evalue_threshold'     => 1e-6,
    'truncated'            => 0,
    'nonsignificant'       => 0,
    'bias_threshold'       => 0.3,
    'gene_overlap_sources' => ['mirbase_gene', 'trnascan_gene'],
  };
}

sub fetch_input {
  my ($self) = @_;
  
  $self->hit_name_counts();
}

sub hit_name_counts {
  my ($self) = @_;
  my $maximum_per_hit_name = $self->param_required('maximum_per_hit_name');
  my $source_logic_name    = $self->param_required('source_logic_name');
  
  my $sql = '
    SELECT
      value as biotype,
      hit_name,
      COUNT(*) as hit_count
    FROM
      dna_align_feature INNER JOIN
      dna_align_feature_attrib USING (dna_align_feature_id) INNER JOIN
      attrib_type USING (attrib_type_id) INNER JOIN
      analysis USING (analysis_id)
    WHERE
      logic_name = ? AND
      code = "rna_gene_biotype"
    GROUP BY
      value,
      hit_name';
  
  my $sth = $self->core_dbh->prepare($sql);
  $sth->execute($source_logic_name);
  
  my %hit_names_to_exclude;
  while (my $results = $sth->fetch) {
    my ($biotype, $hit_name, $hit_count) = @$results;
    
    if (exists $$maximum_per_hit_name{$biotype}) {
      if ($hit_count > $$maximum_per_hit_name{$biotype}) {
        $hit_names_to_exclude{$hit_name} = 1;
      }
    }
  }
  
  $self->param('hit_names_to_exclude', \%hit_names_to_exclude);
}

sub check_thresholds {
  my ($self, $feature) = @_;
  my $evalue_threshold     = $self->param_required('evalue_threshold');
  my $truncated            = $self->param_required('truncated');
  my $nonsignificant       = $self->param_required('nonsignificant');
  my $bias_threshold       = $self->param_required('bias_threshold');
  my $has_structure        = $self->param_required('has_structure');
  my $hit_names_to_exclude = $self->param_required('hit_names_to_exclude');
  
  my $evalue   = $feature->p_value;
  my $hit_name = $feature->hseqname;
  my $attribs  = $feature->get_all_Attributes;
  
  my ($biotype, $cmscan_bias, $cmscan_significant, $cmscan_truncated, $structure);
  foreach my $attrib (@$attribs) {
    if ($attrib->code eq 'rna_gene_biotype') {
      $biotype = $attrib->value;
    }
    if ($attrib->code eq 'cmscan_bias') {
      $cmscan_bias = $attrib->value;
    }
    if ($attrib->code eq 'cmscan_significant') {
      $cmscan_significant = 1;
    }
    if ($attrib->code eq 'cmscan_truncated') {
      $cmscan_truncated = 1;
    }
    if ($attrib->code eq 'ncRNA') {
      $structure = 1;
    }
  }
  
  my $pass_thresholds = 0;
  
  if (defined $evalue) {
    if ($evalue <= $evalue_threshold) {
      $pass_thresholds = 1;
    }
    
    if (exists $$hit_names_to_exclude{$hit_name}) {
      $pass_thresholds = 0;
    }
    
    if ($biotype eq 'misc_RNA') {
      $pass_thresholds = 0;
    }
    
    if (defined $cmscan_bias) {
      if ($cmscan_bias > $bias_threshold) {
        $pass_thresholds = 0;
      }
    }
    
    if (defined $cmscan_truncated) {
      if (! $truncated && $cmscan_truncated) {
        $pass_thresholds = 0;
      }
    }
    
    if (defined $cmscan_significant) {
      if (! $nonsignificant && ! $cmscan_significant) {
        $pass_thresholds = 0;
      }      
    }
    
    if ($has_structure && ! defined $structure) {
      $pass_thresholds = 0;
    }
  }
  
  return $pass_thresholds;
}

1;
