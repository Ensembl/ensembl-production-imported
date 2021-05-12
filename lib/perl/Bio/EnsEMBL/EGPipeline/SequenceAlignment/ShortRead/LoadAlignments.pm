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

package Bio::EnsEMBL::EGPipeline::DNASequenceAlignment::LoadAlignments;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::EGPipeline::Common::Aligner;
use Bio::EnsEMBL::EGPipeline::Common::SAMParser;

sub param_defaults {
  my ($self) = @_;
  
  return {
    'db_type'      => 'otherfeatures',
    'insdc_ids'    => 1,
    'api_load'     => 1,
    'samtools_dir' => undef,
  };
}

sub run {
  my ($self) = @_;
  
  my $db_type     = $self->param_required('db_type');
  my $logic_name  = $self->param_required('logic_name');
  my $insdc_ids   = $self->param_required('insdc_ids');
  my $api_load    = $self->param_required('api_load');
  
  my $dba = $self->get_DBAdaptor($db_type);
  
  my $aa = $dba->get_adaptor('Analysis');
  my $analysis = $aa->fetch_by_logic_name($logic_name);
  
  my $external_db_id = undef;
  if ($insdc_ids) {
    # This is the 'ENA' entry
    $external_db_id = 700;
  }
  
  my $sam_file   = $self->generate_sam();
  my $sam_parser = Bio::EnsEMBL::EGPipeline::Common::SAMParser->new();
  my $alignments = $sam_parser->parse_sam($sam_file);
  
  foreach my $seq_name (keys %$alignments) {
    foreach my $alignment (@{$$alignments{$seq_name}}) {
	    # In Ensembl, 'N' is not supported in the cigar_line, replace with 'I'
	    $alignment->{cigar_line} =~ s/N/I/g;
      
	    # Swap 'D' and 'I' because STAR reports them relative to the
      # hit sequence but we want them relative to the genome sequence.
	    $alignment->{cigar_line} =~ s/I/Z/g;
	    $alignment->{cigar_line} =~ s/[D|N]/I/g;
	    $alignment->{cigar_line} =~ s/Z/D/g;
    }
  }
  
  if ($api_load) {
    $self->api_load($dba, $analysis, $external_db_id, $alignments);
    
  } else {
    my $dbh = $dba->dbc->db_handle();
    my $analysis_id = $analysis->dbID();
    
    my %seq_region_ids;
    my $sa = $dba->get_adaptor('Slice');
    foreach my $slice (@{$sa->fetch_all('toplevel')}) {
      $seq_region_ids{$slice->seq_region_name()} = $slice->dbID();
    }
    
    $self->sql_load($dbh, $analysis_id, $external_db_id, \%seq_region_ids, $alignments);
  }
}

sub api_load {
  my ($self, $dba, $analysis, $external_db_id, $alignments) = @_;
  
  my $sa   = $dba->get_adaptor('Slice');
  my $dafa = $dba->get_adaptor('DnaAlignFeature');
  
  foreach my $seq_name (keys %$alignments) {
    my $slice = $sa->fetch_by_region('toplevel', $seq_name);
    
    if (!defined $slice) {
      $self->warning("Failed to get a slice object for seq_region '$seq_name'\n");
    }
    
    my @dafs;
    foreach my $alignment (@{$$alignments{$seq_name}}) {
      my $daf = Bio::EnsEMBL::DnaDnaAlignFeature->new(
        -slice          => $slice,
        -start          => $alignment->{seq_start},
        -end            => $alignment->{seq_end},
        -strand         => $alignment->{seq_strand},
        -hstart         => $alignment->{hit_start},
        -hend           => $alignment->{hit_end},
        -hstrand        => $alignment->{hit_strand},
        -hseqname       => $alignment->{read_name},
        -cigar_string   => $alignment->{cigar_line},
        -score          => $alignment->{score},
        -percent_id     => $alignment->{perc_id},
        -hcoverage      => $alignment->{hcoverage},
        -analysis       => $analysis,
        -external_db_id => $external_db_id,
      );
      push @dafs, $daf;
    }
    $dafa->store(@dafs);
  }
}

sub sql_load {
  my ($self, $dbh, $analysis_id, $external_db_id, $seq_region_ids, $alignments) = @_;
  
  $dbh->do("LOCK TABLES `dna_align_feature` WRITE;");
  $dbh->do("/*!40000 ALTER TABLE `dna_align_feature` DISABLE KEYS */");
  
  foreach my $seq_name (keys %$alignments) {
    my $seq_region_id = $$seq_region_ids{$seq_name};
    my $insert_sql = '
      INSERT INTO `dna_align_feature` (
        seq_region_id,
        seq_region_start,
        seq_region_end,
        seq_region_strand,
        hit_start,
        hit_end,
        hit_strand,
        hit_name,
        cigar_line,
        analysis_id,
        score,
        perc_ident,
        hcoverage,
        external_db_id
      ) 
      VALUES ';
    my @values_sql;
    
    foreach my $alignment (@{$$alignments{$seq_name}}) {
      my $seq_start  = $alignment->{seq_start};
      my $seq_end    = $alignment->{seq_end};
      my $seq_strand = $alignment->{seq_strand};
      my $hit_start  = $alignment->{hit_start};
      my $hit_end    = $alignment->{hit_end};
      my $hit_strand = $alignment->{hit_strand};
      my $read_name  = $alignment->{read_name};
      my $cigar_line = $alignment->{cigar_line};
      my $score      = $alignment->{score} || 'NULL';
      my $perc_id    = $alignment->{perc_id};
      my $hcoverage  = $alignment->{hcoverage};
      
      push @values_sql,
        "(
          $seq_region_id,
          $seq_start,
          $seq_end,
          $seq_strand,
          $hit_start,
          $hit_end,
          $hit_strand,
          '$read_name',
          '$cigar_line',
          $analysis_id,
          $score,
          $perc_id,
          $hcoverage,
          $external_db_id
        )";
    }
    $dbh->do($insert_sql.join(',', @values_sql));
  }
  
  $dbh->do("/*!40000 ALTER TABLE `dna_align_feature` ENABLE KEYS */");
  $dbh->do("UNLOCK TABLES");
}

sub generate_sam {
  my ($self) = @_;
  
  my $samtools_dir = $self->param('samtools_dir');
  my $bam_file     = $self->param_required('merged_bam_file');
  my $genome_file  = $self->param_required('genome_file');
  
  my $aligner = Bio::EnsEMBL::EGPipeline::Common::Aligner->new(
    -samtools_dir => $samtools_dir,
  );
  my $calmd_file = "$bam_file.calmd";
  my $sam_file = "$bam_file.sam";
  
  my $samtools = 'samtools';
  if (defined $samtools_dir) {
    $samtools = "$samtools_dir/$samtools"
  }
  
  # Can assume that BAM file has already been sorted by MergeBamSet.
  my $cmd = "$samtools calmd -ue $bam_file $genome_file > $calmd_file";
  system($cmd) == 0 || $self->throw("Cannot execute $cmd");

  $cmd = "$samtools view -h $calmd_file > $sam_file";
  system($cmd) == 0 || $self->throw("Cannot execute $cmd");
  
  return $sam_file;
}

1;
