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

Bio::EnsEMBL::EGPipeline::RNAFeatures::FetchAlignments

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::RNAFeatures::FetchAlignments;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub param_defaults {
  my $self = shift @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'db_type' => 'core',
  };
}

sub run {
  my ($self) = @_;
  my $species       = $self->param_required('species');
  my $run_cmscan    = $self->param_required('run_cmscan');
  my $run_trnascan  = $self->param_required('run_trnascan');
  my $load_mirbase  = $self->param_required('load_mirbase');
  my $alignment_dir = $self->param_required('alignment_dir');
  
  my $dba = $self->get_DBAdaptor($self->param_required('db_type'));
  my $dbh = $dba->dbc->db_handle();
  
  if ($run_cmscan) {
    my $rfam_logic_name   = $self->param_required('rfam_logic_name');
    my $cmscan_cm_file    = $self->param_required('cmscan_cm_file');
    my $cmscan_logic_name = $self->param_required('cmscan_logic_name');
    
    my $logic_name;
    if (! exists $$cmscan_cm_file{$species} && ! exists $$cmscan_cm_file{'all'}) {
      $logic_name = $rfam_logic_name;
    } elsif (exists $$cmscan_logic_name{$species}) {
      $logic_name = $$cmscan_logic_name{$species};
    } elsif (exists $$cmscan_logic_name{'all'}) {
      $logic_name = $$cmscan_logic_name{'all'};
    } else {
      $logic_name = 'cmscan_custom';
    }
    
    my $cmscan_file = "$alignment_dir/cmscan.txt";
    open(my $fh, '>', $cmscan_file) or $self->throw("Failed to open $cmscan_file: $!");
    $self->fetch_alignments($dbh, $logic_name, 'evalue', $fh);
    close($fh);
  }
  
  if ($run_trnascan) {    
    my $logic_name = $self->param_required('trnascan_logic_name');
    
    my $trnascan_file = "$alignment_dir/trnascan.txt";
    open(my $fh, '>', $trnascan_file) or $self->throw("Failed to open $trnascan_file: $!");
    $self->fetch_alignments($dbh, $logic_name, 'score', $fh);
    close($fh);
  }
  
  if ($load_mirbase) {    
    my $logic_name = $self->param_required('mirbase_logic_name');
    
    my $mirbase_file = "$alignment_dir/mirbase.txt";
    open(my $fh, '>', $mirbase_file) or $self->throw("Failed to open $mirbase_file: $!");
    $self->fetch_alignments($dbh, $logic_name, 'score', $fh);
    close($fh);
  }
  $dba->dbc->disconnect_if_idle();
}

sub fetch_alignments {
  my ($self, $dbh, $logic_name, $value_column, $fh) = @_;
  
  my $sql = "
    SELECT
      $value_column,
      m.meta_value AS species,
      dafa.value AS biotype,
      daf.hit_name
    FROM meta m, analysis a
    INNER JOIN dna_align_feature daf USING (analysis_id)
    INNER JOIN dna_align_feature_attrib dafa USING (dna_align_feature_id)
    INNER JOIN attrib_type at USING (attrib_type_id)
    WHERE a.logic_name = ?
    AND m.meta_key = 'species.display_name'
    AND at.code = 'rna_gene_biotype'
    ORDER BY species, biotype, hit_name
  ;";
  
  my $sth = $dbh->prepare($sql);
  $sth->execute($logic_name);
  
  while (my @row = $sth->fetchrow_array) {
    print $fh join("\t", @row), "\n";
  }
}

1;
