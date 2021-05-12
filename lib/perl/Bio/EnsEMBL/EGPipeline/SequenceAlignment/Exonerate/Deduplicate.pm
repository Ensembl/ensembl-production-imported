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

package Bio::EnsEMBL::EGPipeline::SequenceAlignment::Exonerate::Deduplicate;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub param_defaults {
  my $self = shift @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'db_type' => 'otherfeatures',
  };
}

sub run {
  my ($self) = @_;
  
  my $dba = $self->get_DBAdaptor($self->param('db_type'));
  my $dbc = $dba->dbc();
  
  my @exon_sql = $self->deduplicate_exon;
  foreach my $exon_sql (@exon_sql) {
    $dbc->do($exon_sql);
  }
  
  my @transcript_sql = $self->deduplicate_transcript;
  foreach my $transcript_sql (@transcript_sql) {
    $dbc->do($transcript_sql);
  }
  
  my @gene_sql = $self->deduplicate_gene;
  foreach my $gene_sql (@gene_sql) {
    $dbc->do($gene_sql);
  }
}

sub deduplicate_exon {
  my ($self) = @_;
  
  my @sql = (
  'create temporary table tmp_singleton_exon as
    select
      min(e.exon_id) as exon_id,
      e.seq_region_id,
      e.seq_region_start,
      e.seq_region_end,
      e.seq_region_strand
    from
      exon e
    group by
      e.seq_region_id,
      e.seq_region_start,
      e.seq_region_end,
      e.seq_region_strand;',
  
  'create unique index exon_id_idx on tmp_singleton_exon (exon_id);',
  'create index seq_region_idx on tmp_singleton_exon (seq_region_id, seq_region_start);',
  
  'create temporary table tmp_duplicate_exon as
    select
      e.exon_id,
      e.seq_region_id,
      e.seq_region_start,
      e.seq_region_end,
      e.seq_region_strand
    from
      exon e left outer join
      tmp_singleton_exon using (exon_id)
    where tmp_singleton_exon.exon_id is null;',
  
  'create index seq_region_idx on tmp_duplicate_exon (seq_region_id, seq_region_start);',
  
  'create temporary table tmp_exon_map as
    select
      tmp_singleton_exon.exon_id,
      tmp_duplicate_exon.exon_id as duplicate_exon_id
    from
      tmp_singleton_exon inner join
      tmp_duplicate_exon using (seq_region_id, seq_region_start, seq_region_end, seq_region_strand);',
  
  'create unique index duplicate_exon_id_idx on tmp_exon_map (duplicate_exon_id);',
  
  'update
    exon_transcript et inner join
    tmp_exon_map em on et.exon_id = em.duplicate_exon_id
  set et.exon_id = em.exon_id;',
  
  'create temporary table tmp_supporting_feature as
  select exon_id, feature_type, feature_id from supporting_feature;',
  
  'update
    tmp_supporting_feature sf inner join
    tmp_exon_map em on sf.exon_id = em.duplicate_exon_id
  set sf.exon_id = em.exon_id;',
  
  'truncate table supporting_feature;',
  
  'insert into supporting_feature
    select distinct exon_id, feature_type, feature_id from tmp_supporting_feature;',
  
  'delete e.* from
    exon e inner join
    tmp_duplicate_exon using (exon_id);',
  
  'drop temporary table tmp_singleton_exon;',
  'drop temporary table tmp_duplicate_exon;',
  'drop temporary table tmp_exon_map;',
  'drop temporary table tmp_supporting_feature;',
  );
  
  return @sql;
}

sub deduplicate_transcript {
  my ($self) = @_;
  
  my @sql = (
  'SET @@group_concat_max_len=100000;',
  
  'create temporary table tmp_transcript as
    select
      t.transcript_id,
      t.seq_region_id,
      t.seq_region_start,
      t.seq_region_end,
      t.seq_region_strand,
      group_concat(et.exon_id SEPARATOR ",") as exon_ids
    from
      transcript t inner join
      exon_transcript et using (transcript_id)
    group by
      t.transcript_id,
      t.seq_region_id,
      t.seq_region_start,
      t.seq_region_end,
      t.seq_region_strand;',
  
  'create temporary table tmp_singleton_transcript as
    select
      min(t.transcript_id) as transcript_id,
      t.seq_region_id,
      t.seq_region_start,
      t.seq_region_end,
      t.seq_region_strand,
      t.exon_ids
    from
      tmp_transcript t
    group by
      t.seq_region_id,
      t.seq_region_start,
      t.seq_region_end,
      t.seq_region_strand,
      t.exon_ids;',
  
  'create unique index transcript_id_idx on tmp_singleton_transcript (transcript_id);',
  'create index seq_region_idx on tmp_singleton_transcript (seq_region_id, seq_region_start);',
  
  'create temporary table tmp_duplicate_transcript as
    select
      t.transcript_id,
      t.seq_region_id,
      t.seq_region_start,
      t.seq_region_end,
      t.seq_region_strand,
      t.exon_ids
    from
      tmp_transcript t left outer join
      tmp_singleton_transcript using (transcript_id)
    where
      tmp_singleton_transcript.transcript_id is null;',
  
  'create index seq_region_idx on tmp_duplicate_transcript (seq_region_id, seq_region_start);',
  
  'create temporary table tmp_transcript_map as
    select
      tmp_singleton_transcript.transcript_id,
      tmp_duplicate_transcript.transcript_id as duplicate_transcript_id
    from
      tmp_singleton_transcript inner join
      tmp_duplicate_transcript using (seq_region_id, seq_region_start, seq_region_end, seq_region_strand, exon_ids);',
  
  'create unique index duplicate_transcript_id_idx on tmp_transcript_map (duplicate_transcript_id);',
  
  'create temporary table tmp_transcript_supporting_feature as
  select transcript_id, feature_type, feature_id from transcript_supporting_feature;',
  
  'update
    tmp_transcript_supporting_feature tsf inner join
    tmp_transcript_map tm on tsf.transcript_id = tm.duplicate_transcript_id
  set tsf.transcript_id = tm.transcript_id;',
  
  'truncate table transcript_supporting_feature;',
  
  'insert into transcript_supporting_feature
    select distinct transcript_id, feature_type, feature_id from tmp_transcript_supporting_feature;',
  
  'delete t.*, et.* from
    transcript t inner join
    exon_transcript et using (transcript_id) inner join
    tmp_duplicate_transcript using (transcript_id);',
  
  'delete g.* from
    gene g left outer join
    transcript t using (gene_id)
  where t.gene_id is null;',
  
  'drop temporary table tmp_transcript;',
  'drop temporary table tmp_singleton_transcript;',
  'drop temporary table tmp_duplicate_transcript;',
  'drop temporary table tmp_transcript_map;',
  'drop temporary table tmp_transcript_supporting_feature;',
  );
  
  return @sql;
}

sub deduplicate_gene {
  my ($self) = @_;
  
  my @sql = (
  'create temporary table tmp_singleton_gene as
    select
      min(g.gene_id) as gene_id,
      g.biotype,
      g.analysis_id,
      g.seq_region_id,
      g.seq_region_start,
      g.seq_region_end,
      g.seq_region_strand
    from
      gene g
    group by
      g.biotype,
      g.analysis_id,
      g.seq_region_id,
      g.seq_region_start,
      g.seq_region_end,
      g.seq_region_strand;',

  'create unique index gene_id_idx on tmp_singleton_gene (gene_id);',
  'create index seq_region_idx on tmp_singleton_gene (seq_region_id, seq_region_start);',

  'create temporary table tmp_duplicate_gene as
    select
      g.gene_id,
      g.biotype,
      g.analysis_id,
      g.seq_region_id,
      g.seq_region_start,
      g.seq_region_end,
      g.seq_region_strand
    from
      gene g left outer join
      tmp_singleton_gene using (gene_id)
    where tmp_singleton_gene.gene_id is null;',

  'create index seq_region_idx on tmp_duplicate_gene (seq_region_id, seq_region_start);',

  'create temporary table tmp_gene_map as
    select
      tmp_singleton_gene.gene_id,
      tmp_duplicate_gene.gene_id as duplicate_gene_id
    from
      tmp_singleton_gene inner join
      tmp_duplicate_gene using (analysis_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand);',

  'create unique index duplicate_gene_id_idx on tmp_gene_map (duplicate_gene_id);',

  'update
    transcript t inner join
    tmp_gene_map tm on t.gene_id = tm.duplicate_gene_id
  set t.gene_id = tm.gene_id;',

  'delete g.* from
    gene g inner join
    tmp_duplicate_gene using (gene_id);',

  'drop temporary table tmp_singleton_gene;',
  'drop temporary table tmp_duplicate_gene;',
  'drop temporary table tmp_gene_map;',
  );
  
  return @sql;
}

1;
