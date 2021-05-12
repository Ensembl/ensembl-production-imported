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

Bio::EnsEMBL::EGPipeline::RNAFeatures::CreateGenes

=cut

package Bio::EnsEMBL::EGPipeline::RNAFeatures::CreateGenes;

use strict;
use warnings;
use feature 'say';

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::DBEntry;

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub param_defaults {
  my ($self) = @_;
  return {
    'gene_source'          => undef,
    'stable_id_type'       => 'eg',
    'allow_repeat_overlap' => 0,
    'allow_coding_overlap' => 0,
    'compare_field'        => 'evalue',
    'gene_overlap_sources' => [],
  };
}

sub run {
  my ($self) = @_;
  my $source_logic_name    = $self->param_required('source_logic_name');
  my $allow_repeat_overlap = $self->param_required('allow_repeat_overlap');
  my $allow_coding_overlap = $self->param_required('allow_coding_overlap');
  my $compare_field        = $self->param_required('compare_field');
  my $gene_overlap_sources = $self->param_required('gene_overlap_sources');
  my $id_db                = $self->param_required('id_db');
  
  my $dba  = $self->core_dba();
  my $dafa = $dba->get_adaptor('DnaAlignFeature');
  my $ra   = $dba->get_adaptor('RepeatFeature');
  my $ta   = $dba->get_adaptor('Transcript');
  my $aa   = $dba->get_adaptor('Analysis');
  my $ga   = $dba->get_adaptor('Gene');
  my $dbea = $dba->get_adaptor('DBEntry');
  my $tsfa = $dba->get_adaptor('TranscriptSupportingFeature');
  my $ida  = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%$id_db);
  
  my $gene_source = $self->generate_source($dba);
  
  my @features = @{ $dafa->fetch_all_by_logic_name($source_logic_name) };
  
  my ($overlap_count, $repeat_count, $coding_exon_count,
      $rna_gene_count, $thresholds_count,
      $update_count, $create_count) = (0, 0, 0, 0, 0, 0, 0);
  
  # Collate genes, so that we can either push all of them or
  # none of them (in the event of an error that causes the module to fail).
  my @new_genes = ();
  my @updated_genes = ();
  
  FEATURE: foreach my $feature (@features) {
    if (! defined $feature->db_name) {
      my $dbid = $feature->dbID;
      $self->throw("Feature (ID=$dbid) lacks an external_db reference.");
    }
    
    if ($self->dna_align_feature_overlap($dafa, $feature, $compare_field)) {
      $overlap_count++;
      next FEATURE;
    }
    
    if (!$allow_repeat_overlap) {
      my ($overlap, undef) = $self->repeat_feature_overlap($ra, $feature);
      if ($overlap) {
        $repeat_count++;
        next FEATURE;
      }
    }
    
    if (!$allow_coding_overlap) {
      my ($overlap, undef) = $self->coding_exon_overlap($ta, $feature);
      if ($overlap) {
        $coding_exon_count++;
        next FEATURE;
      }
    }
    
    if ($self->check_thresholds($feature)) {
      my ($overlap, $within, $rna_transcript) =
        $self->rna_gene_overlap($ta, $feature, $gene_overlap_sources);
   
      if ($overlap) {
        $rna_gene_count++;
      } elsif ($within) {
        push @updated_genes, [$feature, $rna_transcript];
        $update_count++;
      } else {
        push @new_genes, $feature;
        $create_count++;
      }
    } else {
      $thresholds_count++;
    }
  }
  
  foreach my $new_gene (@new_genes) {
    $self->create_gene($dba, $aa, $ga, $dbea, $ida, $gene_source, $new_gene);
  }
  
  foreach my $updated_gene (@updated_genes) {
    $self->update_gene($aa, $ga, $dbea, $tsfa, @$updated_gene);
  }
  
  my $msg = "$create_count $source_logic_name alignments were converted to genes, ".
            "$overlap_count overlapping alignments were ignored";
  
  if ($repeat_count) {
    $msg .= ", $repeat_count alignments overlapped a repeat feature";
  }
  
  if ($coding_exon_count) {
    $msg .= ", $coding_exon_count alignments overlapped a protein-coding region";
  }
  
  if ($rna_gene_count) {
    $msg .= ", $rna_gene_count alignments overlapped an existing RNA gene";
  }
  
  if ($update_count) {
    $msg .= ", $update_count alignments were merged with an existing RNA gene";
  }
  
  if ($thresholds_count) {
    $msg .= ", $thresholds_count alignments failed QC thresholds";
  }
  
  $self->warning($msg);
}

sub generate_source {
  my ($self, $dba) = @_;
  my $gene_source = $self->param('gene_source');
  
  if (! defined $gene_source) {
    $gene_source = $dba->get_MetaContainer->get_division();
    $gene_source =~ s/([a-z])([A-Z])/$1\_$2/;
    $gene_source = 'Ensembl' unless $gene_source;
  }
  
  return $gene_source;
}

sub dna_align_feature_overlap {
  my ($self, $dafa, $feature, $compare) = @_;
  
  # Check overlaps with other dna_align_features; if a feature with the same
  # hit_name exists in the overlapping set, only make a gene if it
  # has the lowest E-value (default) or highest score.
  
  $compare = 'evalue' unless defined $compare;
  my $overlap = 0;
  
  my @dafs = @{ $dafa->fetch_all_nearest_by_Feature(-FEATURE => $feature, -RANGE => 0) };
  foreach my $daf (@dafs) {
    my $dna_align_feature = $$daf[0];
    if ($feature->dbID != $dna_align_feature->dbID) {
      if ($feature->analysis->logic_name eq $dna_align_feature->analysis->logic_name) {
        if ($feature->hseqname eq $dna_align_feature->hseqname) {
          if (
            ($compare eq 'evalue' && $feature->p_value > $dna_align_feature->p_value) ||
            ($compare eq 'score'  && $feature->score < $dna_align_feature->score)
          ) {
              $self->warning("dna_align_feature (ID:".$dna_align_feature->dbID.") overlaps dna_align_feature (ID:".$feature->dbID.") and has a smaller E-value");
              $overlap = 1;
              last;
          }
        }
      }
    }
  }
  
  return $overlap;
}

sub repeat_feature_overlap {
  my ($self, $ra, $feature) = @_;

  my $f_start = $feature->seq_region_start;
  my $f_end   = $feature->seq_region_end;
  
  my ($overlap, $within) = (0, 0);
  my @repeats = @{ $ra->fetch_all_nearest_by_Feature(-FEATURE => $feature, -RANGE => 0) };
  foreach my $repeat (@repeats) {
    my $repeat_feature = $$repeat[0];
    my $rf_start = $repeat_feature->seq_region_start;
    my $rf_end   = $repeat_feature->seq_region_end;

    if (between($f_start, $rf_start, $rf_end) && between($f_end, $rf_start, $rf_end)) {
      $self->warning("Repeat (ID:".$repeat_feature->dbID.") contains dna_align_feature (ID:".$feature->dbID.")");
      $within = 1;
    } else {
      $self->warning("Repeat (ID:".$repeat_feature->dbID.") overlaps dna_align_feature (ID:".$feature->dbID.")");
      $overlap = 1;
    }
  }  
  
  return ($overlap, $within);
}

sub coding_exon_overlap {
  my ($self, $ta, $feature) = @_;
  
  my $f_start = $feature->seq_region_start;
  my $f_end   = $feature->seq_region_end;
  
  my ($overlap, $within) = (0, 0);
  my @transcripts = @{ $ta->fetch_all_nearest_by_Feature(-FEATURE => $feature, -RANGE => 0) };
  foreach my $t (@transcripts) {
    # Need a fresh copy of the transcripts in order for its slice
    # to encompass the whole transcript.
    my $transcript = $ta->fetch_by_stable_id($$t[0]->stable_id);
    
    foreach my $exon (@{ $transcript->get_all_Exons() }) {
      if (defined $exon->coding_region_start($transcript)) {
        my $e_start = $exon->coding_region_start($transcript);
        my $e_end   = $exon->coding_region_end($transcript);

        if (between($f_start, $e_start, $e_end) && between($f_end, $e_start, $e_end)) {
          $self->warning("Exon (ID:".$exon->dbID.") contains dna_align_feature (ID:".$feature->dbID.")");
          $within = 1;
        } elsif (between($f_start, $e_start, $e_end) || between($f_end, $e_start, $e_end) || between($e_start, $f_start, $f_end)) {
          $self->warning("Exon (ID:".$exon->dbID.") overlaps dna_align_feature (ID:".$feature->dbID.")");
          $overlap = 1;
        }
      }
    }
  }  
  
  return ($overlap, $within);
}

sub rna_gene_overlap {
  my ($self, $ta, $feature, $logic_names) = @_;
  
  # Check overlaps with other RNA genes; if the feature is entirely within
  # an existing gene, it is merged, otherwise it is ignored.
  
  my $f_start = $feature->seq_region_start;
  my $f_end   = $feature->seq_region_end;
  
  my %logic_names = map { $_ => 1 } @$logic_names;
  
  my ($overlap, $within, $rna_transcript) = (0, 0, undef);
  my @transcripts = @{ $ta->fetch_all_nearest_by_Feature(-FEATURE => $feature, -RANGE => 0) };
  foreach my $t (@transcripts) {
    my $transcript = $$t[0]->transfer($feature->slice->seq_region_Slice());

    if (exists $logic_names{$transcript->analysis->logic_name}) {
      my $t_start = $transcript->seq_region_start;
      my $t_end   = $transcript->seq_region_end;

      if (between($f_start, $t_start, $t_end) && between($f_end, $t_start, $t_end)) {
        $self->warning("Transcript (ID:".$transcript->dbID.") contains dna_align_feature (ID:".$feature->dbID.")");
        $within = 1;
        $rna_transcript = $transcript;
      } else {
        $self->warning("Transcript (ID:".$transcript->dbID.") overlaps dna_align_feature (ID:".$feature->dbID.")");
        $overlap = 1;
      }
    }
  }  
  
  return ($overlap, $within, $rna_transcript);
}

sub between {
  my ($a, $x, $y) = @_;
  my $between = 0;
  if ($a >= $x && $a <= $y) {
    $between = 1;
  }
  return $between;
}

sub create_gene {
  my ($self, $dba, $aa, $ga, $dbea, $ida, $gene_source, $feature) = @_;
  my $target_logic_name = $self->param_required('target_logic_name');
  my $stable_id_type    = $self->param_required('stable_id_type');
  
  my $analysis = $aa->fetch_by_logic_name($target_logic_name);
  
  my ($gene_stable_id, $transcript_stable_id) =
    $self->generate_stable_id($dba, $ga, $ida, $stable_id_type, $feature);
  
  my $gene = $self->new_gene($feature, $gene_stable_id, $gene_source);
  $gene->analysis($analysis);
  
  my $transcript = $self->new_transcript($feature, $transcript_stable_id, $gene_source);
  $transcript->analysis($analysis);
  
  my $exon = $self->new_exon($feature, $gene_stable_id);
  
  $transcript->add_Exon($exon);
  $gene->add_Transcript($transcript);
  
  $ga->store($gene);
  
  $gene->canonical_transcript($transcript);
  $ga->update($gene);
            
  my $xref = $self->add_xref($dbea, $feature, $gene, $analysis);
  
  $gene->display_xref($xref);
  $ga->update($gene);
}

sub update_gene {
  my ($self, $aa, $ta, $dbea, $tsfa, $feature, $transcript) = @_;
  my $target_logic_name = $self->param_required('target_logic_name');
  
  my $analysis = $aa->fetch_by_logic_name($target_logic_name);
  my $gene     = $transcript->get_Gene();
  
  $tsfa->store($transcript->dbID, [$feature]);
  
  $self->add_xref($dbea, $feature, $gene, $analysis);
}

sub generate_stable_id {
  my ($self, $dba, $ga, $ida, $stable_id_type, $feature) = @_;
  my ($gene_stable_id, $transcript_stable_id);
  
  if ($stable_id_type eq 'eg') {
    my $location_id = $self->feature_location_id($feature);
    
    my $select_sql = "SELECT id FROM ena_identifiers WHERE ena_id = ?";
    my $insert_sql = "INSERT INTO ena_identifiers (ena_id, type) VALUES (?, 'gene')";
    
    my $sth = $ida->dbc->prepare($select_sql);
    $sth->execute($location_id);
    my ($id) = $sth->fetchrow_array();
    
    if (!$id) {    
      $sth = $ida->dbc->prepare($insert_sql);
      $sth->execute($location_id);
      ($id) = $ida->dbc->db_handle->last_insert_id(undef, undef, undef, undef);
    }
    if (length($id) > 9) {
      $self->throw("Numeric portion of identifier exceeds 9 characters");
    }
    
    $gene_stable_id  = 'ENSRNA';
    $gene_stable_id .= substr(($id + 1e9), 1);
    $transcript_stable_id = "$gene_stable_id-T1";
    
  } elsif ($stable_id_type eq 'vb') {
    my $genes = $ga->fetch_all_by_Slice($feature->slice);
    my $exists = 0;
    foreach my $gene (@$genes) {
      if (
        $gene->seq_region_start  == $feature->start &&
        $gene->seq_region_end    == $feature->end &&
        $gene->seq_region_strand == $feature->strand
      ) {
        $gene_stable_id = $gene->stable_id;
        $exists = 1;
        last;
      }
    }
    
    if (!$exists) {
      my $prefix = $dba->get_MetaContainer->single_value_by_key('species.stable_id_prefix');
      my $select_sql = "SELECT max(stable_id) FROM gene WHERE stable_id LIKE '$prefix%'";
      my $sth = $dba->dbc->prepare($select_sql);
      $sth->execute();
      my ($id) = $sth->fetchrow_array();
      
      if (!$id) {
        $gene_stable_id = $prefix.'000001';
      } else {
        $id =~ s/$prefix//;
        my $old_length = length($id);
        $id++;
        my $new_length = length($id);
        my $leading_zeros = ($old_length - $new_length) x '0' || '';
        
        $gene_stable_id = "$prefix$leading_zeros$id";
      }
    }
    
    $transcript_stable_id = "$gene_stable_id-RA";
    
  } else {
    $self->throw("Unknown stable_id_type: $stable_id_type");
  }
  
  return ($gene_stable_id, $transcript_stable_id);
}

sub feature_location_id {
  my ($self, $feature) = @_;
  my $species = $self->param_required('species');
  
  my $accessions = $feature->get_all_Attributes('rfam_accession');
  my $accession  = $feature->hseqname;
  if (scalar(@$accessions) == 1) {
    $accession = $$accessions[0]->value;
  }
  
  my @location_id = (
    $accession,
    $species,
    $feature->coord_system_name,
    $feature->seq_region_name,
    $feature->start,
    $feature->end,
    $feature->strand,
  );
  my $location_id = join(":", @location_id);
  
  return $location_id;
}

sub new_gene {
  my ($self, $feature, $stable_id, $source) = @_;
  
  my $biotypes = $feature->get_all_Attributes('rna_gene_biotype');
  my $biotype  = 'ncRNA';
  if (scalar(@$biotypes) == 1) {
    $biotype = $$biotypes[0]->value;
  }
  
  my $descriptions = $feature->get_all_Attributes('description');
  my $description  = undef;
  if (scalar(@$descriptions) == 1) {
    $description = $$descriptions[0]->value;
  }
  
  my $gene = Bio::EnsEMBL::Gene->new
  (
    -stable_id     => $stable_id,
    -description   => $description,
    -biotype       => $biotype,
    -source        => $source,
    -status        => 'NOVEL',
    -created_date  => time,
    -modified_date => time,
  );
  $gene->version(undef);
  
  return $gene;
}

sub new_transcript {
  my ($self, $feature, $stable_id, $source) = @_;
  
  my $biotypes = $feature->get_all_Attributes('rna_gene_biotype');
  my $biotype  = 'ncRNA';
  if (scalar(@$biotypes) == 1) {
    $biotype = $$biotypes[0]->value;
  }
  
  my $structures = $feature->get_all_Attributes('ncRNA');
  my $structure  = undef;
  if (scalar(@$structures) == 1) {
    $structure = $$structures[0]->value;
  }
  
  my $transcript = Bio::EnsEMBL::Transcript->new
  (
    -stable_id     => $stable_id,
    -start         => $feature->start,
    -end           => $feature->end,
    -strand        => $feature->strand,
    -slice         => $feature->slice,
    -biotype       => $biotype,
    -source        => $source,
    -status        => 'NOVEL',
    -created_date  => time,
    -modified_date => time,
  );
  $transcript->version(undef);
  $transcript->add_supporting_features($feature);
  
  if (defined $structure) {
    my $attrib = Bio::EnsEMBL::Attribute->new(
      -CODE  => 'ncRNA',
      -VALUE => $structure
    );
    $transcript->add_Attributes($attrib);
  }
  
  return $transcript;
}

sub new_exon {
  my ($self, $feature, $gene_stable_id) = @_;
  
  my $exon = Bio::EnsEMBL::Exon->new
  (
    -stable_id       => "$gene_stable_id-E1",
    -start           => $feature->start,
    -end             => $feature->end,
    -strand          => $feature->strand,
    -slice           => $feature->slice,
    -phase           => -1,
    -end_phase       => -1,
    -is_constitutive => 1,
    -created_date    => time,
    -modified_date   => time,
  );
  $exon->version(undef);
  
  return $exon;
}

sub add_xref {
  my ($self, $dbea, $feature, $gene, $analysis) = @_;
  
  my $hit_name = $feature->hseqname;
  my $db_name  = $feature->db_name;
  
  my $accessions = $feature->get_all_Attributes('rfam_accession');
  my $accession  = $hit_name;
  if (scalar(@$accessions) == 1) {
    $accession = $$accessions[0]->value;
  }
  
  my $descriptions = $feature->get_all_Attributes('description');
  my $description  = undef;
  if (scalar(@$descriptions) == 1) {
    $description = $$descriptions[0]->value;
  }
  
  my $xref = Bio::EnsEMBL::DBEntry->new
  (
    -dbname      => $db_name,
    -primary_id  => $accession,
    -display_id  => $hit_name,
    -description => $description,
  );
  $xref->analysis($analysis);
  
  $dbea->store($xref, $gene->dbID, 'Gene');
  
  return $xref;
}

1;
