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

Bio::EnsEMBL::EGPipeline::ProjectGenes::ProjectGenes

=cut

package Bio::EnsEMBL::EGPipeline::ProjectGenes::ProjectGenes;

use strict;
use warnings;
use feature 'say';

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;

use Devel::Cycle;
use Devel::Size qw(size total_size);
$Devel::Size::warn = 0;

use File::Spec::Functions qw(catdir);
use List::Util qw(min max);
use Path::Tiny qw(path);

sub param_defaults {
  my ($self) = @_;

  return {
    logic_name => [],
    filenames  => {
      projected_gff3    => 'projected.gff3',
      projected_cdna    => 'projected_cdna.fa',
      projected_cds     => 'projected_cds.fa',
      projected_pep     => 'projected_pep.fa',
      unprojected_gff3  => 'unprojected.gff3',
      unprojected_cdna  => 'unprojected_cdna.fa',
      unprojected_cds   => 'unprojected_cds.fa',
      unprojected_pep   => 'unprojected_pep.fa',
      exon_report       => 'exons.txt',
      transcript_report => 'transcripts.txt',
      summary_report    => 'summary.txt',
    },
  };
}

sub run {
  my ($self) = @_;
  my $logic_names = $self->param_required('logic_name');
  my $results_dir = $self->param_required('results_dir');
  my $filenames   = $self->param_required('filenames');

  my ($results, $transcript_stats, $summary_stats) =
    $self->initialise_results($results_dir, $filenames);

  $self->project_transcripts($logic_names, $results, $transcript_stats, $summary_stats);
}

sub write_output {
  my ($self) = @_;
  my $results_dir = $self->param_required('results_dir');
  my $filenames   = $self->param_required('filenames');
  
  my $projected_gff3   = catdir($results_dir, $$filenames{'projected_gff3'});
  my $unprojected_gff3 = catdir($results_dir, $$filenames{'unprojected_gff3'});
  
  $self->dataflow_output_id( { 'gff3_file' => $projected_gff3 }, 1);
  $self->dataflow_output_id( { 'gff3_file' => $unprojected_gff3 }, 1);
}

sub fetch_adaptors {
  my ($self) = @_;

  my $from_db_url = $self->param_required('from_db_url');
  my $to_db_url   = $self->param_required('to_db_url');

  # The registry code doesn't like it when you try to
  # instantiate multiple copies of the registry objects,
  # so generate explicit DBA connections instead.

  my $reg = 'Bio::EnsEMBL::Registry';
  my ($from_dba, $to_dba);

  $reg->clear();

  if ($reg->load_registry_from_url("$from_db_url?group=core") == 1) {
    my @dba  = @{ $reg->get_all_DBAdaptors };
    $from_dba = $dba[0]; 
  } else {
    $self->throw("More than one database defined by url $from_db_url");
  }

  my $from_db = {
    -host   => $from_dba->dbc->host,
    -port   => $from_dba->dbc->port,
    -user   => $from_dba->dbc->user,
    -pass   => $from_dba->dbc->pass,
    -dbname => $from_dba->dbc->dbname,
  };

  $reg->clear();

  if ($reg->load_registry_from_url("$to_db_url?group=core") == 1) {
    my @dba  = @{ $reg->get_all_DBAdaptors };
    $to_dba = $dba[0]; 
  } else {
    $self->throw("More than one database defined by url $to_db_url");
  }

  my $to_db = {
    -host   => $to_dba->dbc->host,
    -port   => $to_dba->dbc->port,
    -user   => $to_dba->dbc->user,
    -pass   => $to_dba->dbc->pass,
    -dbname => $to_dba->dbc->dbname,
  };

  $reg->clear();

  $from_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%$from_db);
  $to_dba   = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%$to_db);

  return ($from_dba, $to_dba);
}

sub mapping_target {
  my ($self, $to_dba) = @_;

  my $mca = $to_dba->get_adaptor('MetaContainer');
  my $csa = $to_dba->get_adaptor('CoordSystem');

  my $assembly = $mca->single_value_by_key('assembly.name');
  my @coord_systems = @{$csa->fetch_all_by_version($assembly)};
  @coord_systems = map { $_->name } sort { $a->rank <=> $b->rank } @coord_systems;
  
  return ($assembly, \@coord_systems);
}

sub project_transcripts {
  my ($self, $logic_names, $results, $transcript_stats, $summary_stats) = @_;

  my ($from_dba, $to_dba) = $self->fetch_adaptors();

  my $ta = $from_dba->get_adaptor('Transcript');
  my $sa;

  my ($assembly, $coord_systems) = $self->mapping_target($to_dba);
  
  my $transcripts;
  if (scalar(@$logic_names) == 0) {
    $transcripts = $ta->fetch_all();
  } else {
    foreach my $logic_name (@$logic_names) {
      push @$transcripts, $ta->fetch_all_by_logic_name($logic_name);
    }
  }

  $$summary_stats{'total'} = scalar @$transcripts;

  my %projected_genes;
  my %unprojected_genes;

  my @transcript_ids = map { $_->stable_id } @$transcripts;

  my $counter = 0;

  foreach my $stable_id (@transcript_ids) {
    # There's a memory leak somewhere, but I can't track it down.
    # As a workaround, get fresh adaptors every 100 transcripts.
    if ($counter++ % 100 == 0) {
      say "$counter ".(total_size($ta)/1000000);

      ($from_dba, $to_dba) = $self->fetch_adaptors();
      $ta = $from_dba->get_adaptor('Transcript');
      $sa =   $to_dba->get_adaptor('Slice');
    }

    my $transcript = $ta->fetch_by_stable_id($stable_id);
    my %stats      = map { $_ => 0 } @$transcript_stats;
    $self->transcript_properties($transcript, \%stats);

    my $projected_exons = $self->project_exons($$results{'exon_report'}, $assembly, $coord_systems, $transcript, \%stats);
    my $projected_cds   = $self->project_cds($assembly, $coord_systems, $transcript, \%stats);

    my $coding_transcript;

    if ($stats{'cds_locations'} == 1) {
      ($coding_transcript) = $self->generate_transcript($stable_id, $sa, $projected_cds);
      $stats{'cds_translates'} = $self->check_translation($coding_transcript);
      if (! $stats{'cds_translates'} && $stats{'cds'} == $stats{'cds_complete'}) {
        $stats{'cds_translates'} = $self->shuffle($coding_transcript, $projected_exons, $projected_cds);
        if ($stats{'cds_translates'}) {
          ($coding_transcript) = $self->generate_transcript($stable_id, $sa, $projected_cds);
          $self->warning("Shuffled $stable_id to get translation.");
        }
      }

      $stats{'cdna_edit_dist'}   = $self->levenshtein_distance($transcript->seq->seq, $coding_transcript->seq->seq);
      $stats{'coding_edit_dist'} = $self->levenshtein_distance($transcript->translateable_seq, $coding_transcript->translateable_seq);
    }

    if ($stats{'exons'} == $stats{'exon_complete'}) {
      if ($stats{'exon_locations'} == 1) {
        if ($stats{'biotype'} ne 'protein_coding' || $stats{'cds_translates'}) {
          $stats{'status'} = 'Projected';
          $self->write_gff3_transcript($$results{'projected_gff3'}, $transcript, 0, $projected_exons, $projected_cds, \%projected_genes);
          $self->write_fasta_sequence($$results{'projected_cdna'}, $$results{'projected_cds'}, $$results{'projected_pep'}, $transcript, $stats{'status'});
          $$summary_stats{'transcript_complete'}++;

        } else {
          $stats{'status'} = 'Unprojected (CDS has internal stop codons)';
          $self->write_gff3_transcript($$results{'unprojected_gff3'}, $transcript, 0, $projected_exons, $projected_cds, \%unprojected_genes);
          $self->write_fasta_sequence($$results{'unprojected_cdna'}, $$results{'unprojected_cds'}, $$results{'unprojected_pep'}, $transcript, $stats{'status'});
          $$summary_stats{'transcript_noncoding'}++;

        }

      } elsif ($stats{'cds'} > 0 && $stats{'cds_locations'} == 1) {
        if ($stats{'cds_translates'}) {
          $stats{'status'} = 'Projected CDS (UTR on multiple scaffolds/strands)';
          $self->write_gff3_transcript($$results{'projected_gff3'}, $transcript, 0, $projected_cds, $projected_cds, \%projected_genes);
          $self->write_fasta_sequence($$results{'projected_cdna'}, $$results{'projected_cds'}, $$results{'projected_pep'}, $transcript, $stats{'status'});
          $$summary_stats{'transcript_coding'}++;

        } else {
          $stats{'status'} = 'Unprojected (CDS has internal stop codons)';
          $self->write_gff3_transcript($$results{'unprojected_gff3'}, $transcript, 0, $projected_cds, $projected_cds, \%unprojected_genes);
          $self->write_fasta_sequence($$results{'unprojected_cdna'}, $$results{'unprojected_cds'}, $$results{'unprojected_pep'}, $transcript, $stats{'status'});
          $$summary_stats{'transcript_noncoding'}++;

        }

      } else {
        $stats{'status'} = 'Unprojected (CDS on multiple scaffolds/strands)';
        $self->write_gff3_transcript($$results{'unprojected_gff3'}, $transcript, 1, $projected_exons, $projected_cds, \%unprojected_genes);
        $self->write_fasta_sequence($$results{'unprojected_cdna'}, $$results{'unprojected_cds'}, $$results{'unprojected_pep'}, $transcript, $stats{'status'});
        $$summary_stats{'transcript_split'}++;

      }

    } elsif ($stats{'cds'} > 0 && $stats{'cds'} == $stats{'cds_complete'}) {
      if ($stats{'cds_locations'} == 1) {
        if ($stats{'cds_translates'}) {
          $stats{'status'} = 'Projected CDS (UTR truncated)';
          $self->write_gff3_transcript($$results{'projected_gff3'}, $transcript, 0, $projected_cds, $projected_cds, \%projected_genes);
          $self->write_fasta_sequence($$results{'projected_cdna'}, $$results{'projected_cds'}, $$results{'projected_pep'}, $transcript, $stats{'status'});
          $$summary_stats{'transcript_coding'}++;

        } else {
          $stats{'status'} = 'Unprojected (CDS has internal stop codons)';
          $self->write_gff3_transcript($$results{'unprojected_gff3'}, $transcript, 0, $projected_cds, $projected_cds, \%unprojected_genes);
          $self->write_fasta_sequence($$results{'unprojected_cdna'}, $$results{'unprojected_cds'}, $$results{'unprojected_pep'}, $transcript, $stats{'status'});
          $$summary_stats{'transcript_noncoding'}++;

        }

      } else {
        $stats{'status'} = 'Unprojected (CDS on multiple scaffolds/strands)';
        $self->write_gff3_transcript($$results{'unprojected_gff3'}, $transcript, 1, $projected_exons, $projected_cds, \%unprojected_genes);
        $self->write_fasta_sequence($$results{'unprojected_cdna'}, $$results{'unprojected_cds'}, $$results{'unprojected_pep'}, $transcript, $stats{'status'});
        $$summary_stats{'transcript_split'}++;

      }

    } else {
      if (($stats{'exon_partial'} + $stats{'exon_complete'}) > 0) {
        my $split = $stats{'exon_locations'} > 1 ? 1 : 0;
        $stats{'status'} = 'Unprojected (CDS truncated)';
        $self->write_gff3_transcript($$results{'unprojected_gff3'}, $transcript, $split, $projected_exons, $projected_cds, \%unprojected_genes);
        $self->write_fasta_sequence($$results{'unprojected_cdna'}, $$results{'unprojected_cds'}, $$results{'unprojected_pep'}, $transcript, $stats{'status'});
        $$summary_stats{'transcript_partial'}++;

      } else {
        $stats{'status'} = 'Unprojected';
        $self->write_fasta_sequence($$results{'unprojected_cdna'}, $$results{'unprojected_cds'}, $$results{'unprojected_pep'}, $transcript, $stats{'status'});
        $$summary_stats{'transcript_missing'}++;

      }
    }

    $self->write_transcript_report($$results{'transcript_report'}, $transcript_stats, \%stats);
  }

  $self->write_gff3_genes($$results{'projected_gff3'}, \%projected_genes);
  $self->write_gff3_genes($$results{'unprojected_gff3'}, \%unprojected_genes);

  $self->write_summary_report($$results{'summary_report'}, $summary_stats);

}

sub transcript_properties {
  my ($self, $transcript, $stats) = @_;

  $$stats{'transcript'}  = $transcript->stable_id;
  $$stats{'biotype'}     = $transcript->biotype;
  $$stats{'logic_name'}  = $transcript->analysis->logic_name;
  $$stats{'cdna_length'} = $transcript->length;

  my $coding_length;
  if (defined $transcript->translation) {
    $coding_length = length($transcript->translateable_seq);
  } else {
    $coding_length = $transcript->length;
  }
  $$stats{'coding_length'} = $coding_length;

  if ($transcript->biotype eq 'protein_coding') {
    my $pfs = $transcript->translation->get_all_ProteinFeatures();
    foreach my $pf (@$pfs) {
      $$stats{'protein_features'}++ if $pf->analysis->logic_name ne 'seg';
    }
  }
}

sub project_exons {
  my ($self, $exon_report, $assembly, $coord_systems, $transcript, $stats) = @_;

  my $exons = $transcript->get_all_Exons();
  $$stats{'exons'} = scalar @$exons;

  my %location;

  foreach my $exon (@$exons) {
    my $projection;
    foreach my $coord_system (@$coord_systems) {
      $projection = $exon->project($coord_system, $assembly);
      last if scalar(@$projection) > 0;
    }
    $self->exon_location($stats, \%location, $exon, $projection);
    $self->write_exon_report($exon_report, $exon, $transcript, $projection);
  }

  foreach my $name (keys %location) {
    foreach my $strand (keys %{$location{$name}}) {
      $$stats{'exon_locations'}++;
    }
  }

  return \%location;
}

sub exon_location {
  my ($self, $stats, $location, $exon, $projection) = @_;

  my $seg_length = 0;
  foreach my $seg (@$projection) {
    $seg_length += $seg->from_end - $seg->from_start + 1;

    my $to_slice = $seg->to_Slice();
    my $name = $to_slice->seq_region_name;
    my $strand = $to_slice->strand < 0 ? '-' : '+';
    my $start = $to_slice->start;
    my $end = $to_slice->end;

    push @{$$location{$name}{$strand}}, [$start, $end];
  }

  if ($exon->length == $seg_length) {
    $$stats{'exon_complete'}++;
  } elsif ($seg_length == 0) {
    $$stats{'exon_missing'}++;
  } else {
    $$stats{'exon_partial'}++;
  }
}

sub project_cds {
  my ($self, $assembly, $coord_systems, $transcript, $stats) = @_;

  my $cds_features = $transcript->get_all_CDS();
  $$stats{'cds'} = scalar @$cds_features;

  my %location;

  foreach my $cds (@$cds_features) {
    $cds->slice($transcript->slice);
    my $projection;
    foreach my $coord_system (@$coord_systems) {
      $projection = $cds->project($coord_system, $assembly);
      last if scalar(@$projection) > 0;
    }
    $self->cds_location($stats, \%location, $cds, $projection);
  }

  foreach my $name (keys %location) {
    foreach my $strand (keys %{$location{$name}}) {
      $$stats{'cds_locations'}++;
    }
  }

  return \%location;
}

sub cds_location {
  my ($self, $stats, $location, $cds, $projection) = @_;

  my $phase = $cds->phase;

  my $seg_length = 0;
  foreach my $seg (@$projection) {
    $seg_length += $seg->from_end - $seg->from_start + 1;

    my $to_slice = $seg->to_Slice();
    my $name = $to_slice->seq_region_name;
    my $strand = $to_slice->strand < 0 ? '-' : '+';
    my $start = $to_slice->start;
    my $end = $to_slice->end;

    push @{$$location{$name}{$strand}}, [$start, $end, $phase];
  }

  my $cds_length = $cds->end - $cds->start + 1;
  if ($cds_length == $seg_length) {
    $$stats{'cds_complete'}++;
  } elsif ($seg_length == 0) {
    $$stats{'cds_missing'}++;
  } else {
    $$stats{'cds_partial'}++;
  }
}

sub generate_transcript {
  my ($self, $stable_id, $sa, $projected_cds) = @_;

  my @transcripts;

  foreach my $name (sort keys %{$projected_cds}) {
    foreach my $strand (sort keys %{$$projected_cds{$name}}) {
      my $slice = $sa->fetch_by_region('toplevel', $name);
      my $num_strand = $strand eq '-' ? -1 : 1;

      my @exons;
      foreach my $pos (sort {$$a[0] <=> $$b[0]} @{$$projected_cds{$name}{$strand}}) {
        my $start = $$pos[0];
        my $end = $$pos[1];
        my $phase = scalar(@$pos) == 3 ? $$pos[2] : 0;

        my $exon = Bio::EnsEMBL::Exon->new(
          -slice  => $slice,
          -strand => $num_strand,
          -start  => $start,
          -end    => $end,
          -phase  => $phase =~ tr/12/21/,
        );

        push @exons, $exon;
      }

      if ($strand eq '-') {
        @exons = reverse @exons;
      }

      my $transcript = Bio::EnsEMBL::Transcript->new(
        -stable_id => $stable_id,
        -exons     => \@exons,
      );

      my $translation = Bio::EnsEMBL::Translation->new(
        -start_exon => $exons[0],
        -end_exon   => $exons[-1],
        -seq_start  => 1,
        -seq_end    => $exons[-1]->length,
      );
      $transcript->translation($translation);

      push @transcripts, $transcript;
    }
  }

  return @transcripts;
}

sub check_translation {
  my ($self, $transcript) = @_;

  if (defined $transcript && defined $transcript->translate) {
    my $seq = $transcript->translate->seq;
    if ($seq =~ /\*/) {
      return 0;
    } else {
      return 1;
    }
  } else {
    return 0;
  }
}

sub shuffle {
  my ($self, $new_transcript, $projected_exons, $projected_cds) = @_;

  my $translates = 0;
  my @offsets = (-1, 1, -2, 2);

  foreach my $offset (@offsets) {
    if ($self->shuffle_by_offset($new_transcript, $offset)) {
      $translates = 1;

      $self->shift_locations($projected_exons, $offset);
      $self->shift_locations($projected_cds, $offset);

      last;
    }
  }

  return $translates;
}

sub shuffle_by_offset {
  my ($self, $transcript, $offset) = @_;

  $self->shift_transcript($transcript, $offset);
  if ($transcript->start < 1 || $transcript->end > $transcript->slice->end) {
    $self->shift_transcript($transcript, -1 * $offset);
    return 0;
  }

  if ($self->check_translation($transcript)) {
    return 1;
  } else {
    $self->shift_transcript($transcript, -1 * $offset);
    return 0;
  }
}

sub shift_transcript {
  my ($self, $transcript, $start_adjustment, $end_adjustment) = @_;
  $start_adjustment = -1 unless defined $start_adjustment;
  $end_adjustment = $start_adjustment unless defined $end_adjustment;

  $transcript->start($transcript->start + $start_adjustment);
  $transcript->end($transcript->end + $end_adjustment);

  foreach my $exon (@{$transcript->get_all_Exons}) {
    $exon->start($exon->start + $start_adjustment);
    $exon->end($exon->end + $end_adjustment);
  }
}

sub shift_locations {
  my ($self, $projected_locations, $offset) = @_;

  foreach my $name (sort keys %{$projected_locations}) {
    foreach my $strand (sort keys %{$$projected_locations{$name}}) {
      foreach my $pos (sort {$$a[0] <=> $$b[0]} @{$$projected_locations{$name}{$strand}}) {
        $$pos[0] += $offset;
        $$pos[1] += $offset;
      }
    }
  }
}

sub initialise_results {
  my ($self, $results_dir, $filenames) = @_;

  # These stats are calculated for every exon.
  my @exons_stats = qw(
    transcript exon from_scaffold from_start from_end from_strand
    segments segment_start segment_end to_scaffold to_start to_end to_strand
  );

  # These stats are calculated for every transcript.
  my @transcript_stats = qw(
    transcript status
    logic_name biotype cdna_length coding_length exons cds protein_features
    exon_complete exon_partial exon_missing exon_locations
    cds_complete cds_partial cds_missing cds_locations cds_translates
    cdna_edit_dist coding_edit_dist
  );

  # These stats are calculated across all transcripts.
  my @summary_stats = qw(
    transcript_complete transcript_coding
    transcript_noncoding transcript_partial transcript_split transcript_missing
  );
  my %summary_stats = map { $_ => 0 } @summary_stats;

  # We generate lots of data, so rather than accumulate in variables,
  # write direct to file.

  my %results;
  foreach my $filename (keys %$filenames) {
    path($results_dir)->mkpath;

    my $file = path(catdir($results_dir, $$filenames{$filename}));
    $file->touch;

    if ($filename eq 'exon_report') {
      $file->spew(join("\t", @exons_stats)."\n");
    } elsif ($filename eq 'transcript_report') {
      $file->spew(join("\t", @transcript_stats)."\n");
    } elsif ($filename =~ /(projected_gff3|unprojected_gff3)/) {
      $file->spew("##gff-version 3\n");
    } else {
      $file->spew('');
    }

    $results{$filename} = $file;
  }

  return (\%results, \@transcript_stats, \%summary_stats);
}

sub write_gff3_transcript {
  my ($self, $gff3, $transcript, $split, $projected_exons, $projected_cds, $projected_genes) = @_;

  my $suffix = 1;
  my $gene_id = $transcript->get_Gene->stable_id;
  my $transcript_id = $transcript->stable_id;
  my $biotype = $transcript->biotype;
  my $logic_name = $transcript->analysis->logic_name;

  foreach my $name (sort keys %{$projected_exons}) {
    foreach my $strand (sort keys %{$$projected_exons{$name}}) {
      my $g_id = $gene_id;
      my $t_id = $transcript_id;
      if ($split) {
        $g_id = "$gene_id-$suffix";
        ($t_id = $transcript_id) =~ s/$gene_id/$g_id/;
        $suffix++;
      }

      my @starts;
      my @ends;
      foreach my $pos (@{$$projected_exons{$name}{$strand}}) {
        my ($start, $end) = @$pos;
        push @starts, $start;
        push @ends, $end;

        $gff3->append(join("\t",
          (
            $name,
            $logic_name,
            'exon',
            $start,
            $end,
            '.',
            $strand,
            '.',
            'Parent='.$t_id,
          )
        )."\n");
      }

      if (exists $$projected_cds{$name}{$strand}) {
        foreach my $pos (@{$$projected_cds{$name}{$strand}}) {
          my $start = $$pos[0];
          my $end = $$pos[1];
          my $phase = scalar(@$pos) == 3 ? $$pos[2] : '.';

          $gff3->append(join("\t",
            (
              $name,
              $logic_name,
              'CDS',
              $start,
              $end,
              '.',
              $strand,
              $phase,
              'Parent='.$t_id,
            )
          )."\n");
        }
      }

      # It's possibly (but indicative of a terrible gene model) for two
      # transcripts to map, in their entirety, to different scaffolds...
      foreach my $pos (@{$$projected_genes{$g_id}}) {
        my $gene_name = $$pos[2];
        my $gene_strand =$$pos[3];
        if ($name ne $gene_name) {
          $g_id = "$g_id-multi_scaf";
          last;
        }
        if ($strand ne $gene_strand) {
          $g_id = "$g_id-multi_strand";
          last;
        }
      }

      my $transcript_type = $biotype;
      $transcript_type = 'mRNA' if $biotype eq 'protein_coding';
      $transcript_type = 'transcript' if $biotype eq 'pseudogene';
      my ($t_start, $t_end) = (min(@starts), max(@ends));
      $gff3->append(join("\t",
        (
          $name,
          $logic_name,
          $transcript_type,
          $t_start,
          $t_end,
          '.',
          $strand,
          '.',
          'ID='.$t_id.';Parent='.$g_id,
        )
      )."\n");

      push @{$$projected_genes{$g_id}}, [$logic_name, $biotype, $name, $strand, $t_start, $t_end];
    }
  }
}

sub write_gff3_genes {
  my ($self, $gff3, $projected_genes) = @_;

  foreach my $g_id (sort keys %$projected_genes) {
    my $logic_name;
    my $biotype;
    my $name;
    my $strand;
    my @starts;
    my @ends;
    foreach my $pos (@{$$projected_genes{$g_id}}) {
      $logic_name = $$pos[0];
      $biotype = $$pos[1];
      $name = $$pos[2];
      $strand =$$pos[3];
      push @starts, $$pos[4];
      push @ends, $$pos[5];
    }

    $biotype = 'gene' if $biotype ne 'pseudogene';
    my ($g_start, $g_end) = (min(@starts), max(@ends));
    $gff3->append(join("\t",
      (
        $name,
        $logic_name,
        $biotype,
        $g_start,
        $g_end,
        '.',
        $strand,
        '.',
        'ID='.$g_id,
      )
    )."\n");
  }
}

sub write_fasta_sequence {
  my ($self, $cdna_data, $cds_data, $pep_data, $transcript, $status) = @_;

  my $header = ">".$transcript->stable_id." $status";

  (my $cdna = $transcript->seq()->seq()) =~ s/([\w\-\.]{1,80})/$1\n/gms;
	$cdna_data->append("$header\n$cdna\n");

  (my $cds = $transcript->translateable_seq()) =~ s/([\w\-\.]{1,80})/$1\n/gms;
	$cds_data->append("$header\n$cds\n");

  if (defined $transcript->translation) {
    (my $pep = $transcript->translate()->seq()) =~ s/([\w\-\.]{1,80})/$1\n/gms;
  	$pep_data->append("$header\n$pep\n");
  }
}

sub write_exon_report {
  my ($self, $exon_report, $exon, $transcript, $projection) = @_;

  my $segments = scalar @$projection;
  if ($segments) {
    foreach my $seg (@$projection) {
      my $to_slice = $seg->to_Slice();
      $exon_report->append(join("\t",
        (
          $transcript->stable_id,
          $exon->stable_id,
          $exon->slice->seq_region_name,
          $exon->start,
          $exon->end,
          $exon->strand,
          $segments,
          $seg->from_start,
          $seg->from_end,
          $to_slice->seq_region_name,
          $to_slice->start,
          $to_slice->end,
          $to_slice->strand,
        )
      )."\n");
    }
  } else {
    $exon_report->append(join("\t",
      (
        $transcript->stable_id,
        $exon->stable_id,
        $exon->slice->seq_region_name,
        $exon->start,
        $exon->end,
        $exon->strand,
        $segments,
        '',
        '',
        '',
        '',
        '',
        '',
      )
    )."\n");
  }
}

sub write_transcript_report {
  my ($self, $transcript_report, $stat_names, $stats) = @_;

  $transcript_report->append(join("\t", (map { $$stats{$_} } @$stat_names))."\n");
}

sub write_summary_report {
  my ($self, $summary_report, $summary_stats) = @_;

  my $total       = $$summary_stats{'total'};
  my $projected   = $$summary_stats{'transcript_complete'} +
                    $$summary_stats{'transcript_coding'};
  my $unprojected = $$summary_stats{'transcript_noncoding'} +
                    $$summary_stats{'transcript_partial'} +
                    $$summary_stats{'transcript_split'} +
                    $$summary_stats{'transcript_missing'};

  my @summary = (
    "Transcript Summary",

    "Total transcripts: $total",

    "Projected: ".
      stat_text($projected, $total, 'transcripts'),

    "  With complete coding+utr sequence: ".
      stat_text($$summary_stats{'transcript_complete'}, $projected, 'projected'),

    "  With complete coding sequence: ".
      stat_text($$summary_stats{'transcript_coding'}, $projected, 'projected'),

    "Unprojected: ".
      stat_text($unprojected, $total, 'transcripts'),

    "  With non-translating sequence: ".
      stat_text($$summary_stats{'transcript_noncoding'}, $unprojected, 'unprojected'),

    "  With partial coding sequence: ".
      stat_text($$summary_stats{'transcript_partial'}, $unprojected, 'unprojected'),

    "  With multiple scaffolds/strands: ".
      stat_text($$summary_stats{'transcript_split'}, $unprojected, 'unprojected'),

    "  With no coding sequence: ".
      stat_text($$summary_stats{'transcript_missing'}, $unprojected, 'unprojected'),

  );

  $summary_report->append(join("\n", @summary)."\n");
}

sub stat_text {
  my ($value, $total, $total_type) = @_;
  my $summary;
  if ($total == 0) {
    $summary = "0 (0%";
  } else {
    $summary = "$value (".sprintf("%.1f", $value/$total*100)."% ";
  }
  $summary .= "of $total_type)";
  return $summary;
}

# Algorithm from wikipedia pseudocode.
sub levenshtein_distance {
  my ($self, $s, $t) = @_;
  my ($s_length, $t_length) = (length($s), length($t));

	# degenerate cases
	return (0, $s_length, $t_length) if ($s eq $t);
	return ($t_length, $s_length, $t_length) if ($s_length == 0);
	return ($s_length, $s_length, $t_length) if ($t_length == 0);

  # enable easy access to string characters
  my @s = split(//, $s);
  my @t = split(//, $t);

	# create two work vectors of integer distances
	my (@v0, @v1);

	# initialize v0 (the previous row of distances)
	# this row is A[0][i]: edit distance for an empty s
	# the distance is just the number of characters to delete from t
	for (my $i = 0; $i <= $t_length; $i++) {
		$v0[$i] = $i;
  }

	for (my $i = 0; $i < $s_length; $i++) {
		# calculate v1 (current row distances) from the previous row v0

		# first element of v1 is A[i+1][0]
		#   edit distance is delete (i+1) chars from s to match empty t
		$v1[0] = $i + 1;

		# use formula to fill in the rest of the row
		for (my $j = 0; $j < $t_length; $j++) {
			my $cost = ($s[$i] eq $t[$j]) ? 0 : 1;
			$v1[$j + 1] = min($v1[$j] + 1, $v0[$j + 1] + 1, $v0[$j] + $cost);
		}

		# copy v1 (current row) to v0 (previous row) for next iteration
		for (my $j = 0; $j < scalar(@v0); $j++) {
			$v0[$j] = $v1[$j];
    }
	}

	return $v1[$t_length];
}

1;
