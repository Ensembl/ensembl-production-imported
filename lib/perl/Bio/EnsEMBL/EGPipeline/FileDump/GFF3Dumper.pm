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

package Bio::EnsEMBL::EGPipeline::FileDump::GFF3Dumper;

use strict;
use warnings;
no  warnings 'redefine';
use feature 'say';

use base ('Bio::EnsEMBL::EGPipeline::FileDump::BaseDumper');

use Bio::EnsEMBL::Utils::IO::GFFSerializer;

use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::RepeatFeature;

use Path::Tiny qw(path);

sub param_defaults {
  my ($self) = @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'feature_type'          => ['Gene', 'Transcript'],
    'data_type'             => 'basefeatures',
    'file_type'             => 'gff3',
    'per_chromosome'        => 0,
    'include_scaffold'      => 1,
    'logic_name'            => [],
    'remove_id_prefix'      => 0,
    'relabel_transcript'    => 0,
    'join_align_feature'    => 0,
    'remove_separators'     => 0,
    'dump_level'            => 'toplevel',
    'dump_cs_version'       => undef,
    'include_non_reference' => 0,
  };
}

sub run {
  my ($self) = @_;
  my $species               = $self->param_required('species');
  my $db_type               = $self->param_required('db_type');
  my $out_file              = $self->param_required('out_file');
  my $feature_types         = $self->param_required('feature_type');
  my $per_chromosome        = $self->param_required('per_chromosome');
  my $include_scaffold      = $self->param_required('include_scaffold');
  my $logic_names           = $self->param_required('logic_name');
  my $remove_id_prefix      = $self->param_required('remove_id_prefix');
  my $relabel_transcript    = $self->param_required('relabel_transcript');
  my $join_align_feature    = $self->param_required('join_align_feature');
  my $remove_separators     = $self->param_required('remove_separators');
  my $dump_level            = $self->param_required('dump_level');
  my $dump_cs_version       = $self->param('dump_cs_version');
  my $include_non_reference = $self->param('include_non_reference'); 
  
  my $reg = 'Bio::EnsEMBL::Registry';
  
  my $sa = $reg->get_adaptor($species, $db_type, 'Slice');
  my $slices = $sa->fetch_all($dump_level, $dump_cs_version, $include_non_reference);
  
  my $oa = $reg->get_adaptor('multi', 'ontology', 'OntologyTerm');
  die("Could not load the ontology adaptor") if not $oa;
  
  open(my $out_fh, '>', $out_file) or $self->throw("Cannot open file $out_file: $!");
  my $serializer = Bio::EnsEMBL::Utils::IO::GFFSerializer->new($oa, $out_fh);
  $serializer->print_main_header($slices);
  
  my $mc = $self->core_dba->get_MetaContainer();
  my $provider = $mc->single_value_by_key('provider.name');
  my $assembly = $mc->single_value_by_key('assembly.default');
  $serializer->print_metadata("#genome-build $provider $assembly");
  
  my %adaptors;
  foreach my $feature_type (@$feature_types) {
    $adaptors{$feature_type} = $reg->get_adaptor($species, $db_type, $feature_type);
    if ($feature_type eq 'Transcript') {
      $adaptors{'Exon'} = $reg->get_adaptor($species, $db_type, 'Exon');
    }
  }
  
  my %chr;
  my $has_chromosomes = $self->has_chromosomes($sa);
  
  foreach my $slice (@$slices) {
    if ($include_scaffold) {
      $slice->source($provider) if $provider;
      $serializer->print_feature($slice);
    }
    
    foreach my $feature_type (@$feature_types) {
      my $features = $self->fetch_features($feature_type, $adaptors{$feature_type}, $logic_names, $slice);
      $serializer->print_feature_list($features);
      
      if ($per_chromosome && $has_chromosomes) {
        my $chr_serializer = $self->chr_serializer($oa, $slice, \%chr);
        $chr_serializer->print_feature_list($features);
      }
    }
  }
  
  close($out_fh);
  my $out_files = $self->param('out_files');
  
  foreach my $slice_name (keys %chr) {
    close($chr{$slice_name}{'fh'});
    push @$out_files, $chr{$slice_name}{'file'};
  }
  
  foreach my $out_file (@$out_files) {
    $self->remove_id_prefix($out_file) if $remove_id_prefix;
    $self->relabel_transcript($out_file) if $relabel_transcript;
    $self->join_align_feature($out_file) if $join_align_feature;
    $self->remove_separators($out_file) if $remove_separators;
  }
  
  $self->param('out_files', $out_files);
}

sub fetch_features {
  my ($self, $feature_type, $adaptor, $logic_names, $slice) = @_;
  
  my @features;
  if (scalar(@$logic_names) == 0) {
    @features = @{$adaptor->fetch_all_by_Slice($slice)};
  } else {
    foreach my $logic_name (@$logic_names) {
      my $features;
      if ($feature_type eq 'Transcript') {
        $features = $adaptor->fetch_all_by_Slice($slice, 0, $logic_name);
      } else {
        $features = $adaptor->fetch_all_by_Slice($slice, $logic_name);
      }
      push @features, @$features;
    }
  }
  
  if ($feature_type eq 'Transcript') {
    my $exon_features = $self->exon_features(\@features);
    push @features, @$exon_features;
  }
      
  return \@features;
}

sub exon_features {
  my ($self, $transcripts) = @_;
  
  my @cds_features;
  my @exon_features;
  my @utr_features;
  
  foreach my $transcript (@$transcripts) {
    push @cds_features,  @{ $transcript->get_all_CDS() };
    push @exon_features, @{ $transcript->get_all_ExonTranscripts() };
    
    my $seq_edits = $transcript->get_all_SeqEdits();
    if (! defined $seq_edits || scalar(@$seq_edits) == 0) {
      push @utr_features,  @{ $transcript->get_all_five_prime_UTRs() };
      push @utr_features,  @{ $transcript->get_all_three_prime_UTRs() };
    }
  }
  
  return [@exon_features, @cds_features, @utr_features];
}

sub chr_serializer {
  my ($self, $oa, $slice, $chr) = @_;
  
  my $out_file = $self->param_required('out_file');
  
  my $slice_name;
  if ($slice->karyotype_rank > 0) {
    $slice_name = 'chromosome.'.$slice->seq_region_name;
  } else {
    $slice_name = 'nonchromosomal';
  }
    
  unless (exists $$chr{$slice_name}) {
    (my $chr_file = $out_file) =~ s/([^\.]+)$/$slice_name.$1/;
    open(my $chr_fh, '>', $chr_file) or $self->throw("Cannot open file $chr_file: $!");
    
    my $chr_serializer = Bio::EnsEMBL::Utils::IO::GFFSerializer->new($oa, $chr_fh);
    $chr_serializer->print_main_header([$slice]);
    
    $$chr{$slice_name}{'fh'} = $chr_fh;
    $$chr{$slice_name}{'file'} = $chr_file;
    $$chr{$slice_name}{'serializer'} = $chr_serializer;
  }
  
  return $$chr{$slice_name}{'serializer'};
}

sub remove_id_prefix {
  my ($self, $out_file) = @_;
  
  my $file = path($out_file);
  
  my $data = $file->slurp;
  $data =~ s/(ID|Parent)=\w+:/$1=/gm;
  $file->spew($data);
}

sub relabel_transcript {
  my ($self, $out_file) = @_;
  
  my $file = path($out_file);
  
  my $data = $file->slurp;
  $data =~ s/\ttranscript\t(.*biotype=protein_coding;)/\tmRNA\t$1/gm;
  $data =~ s/\ttranscript\t(.*biotype=(\w*RNA\w*);)/\t$2\t$1/gm;
  $file->spew($data);
}

sub remove_separators {
  my ($self, $out_file) = @_;
  
  my $file = path($out_file);
  
  my $data = $file->slurp;
  $data =~ s/^###\n//m;
  $file->spew($data);
}

sub join_align_feature {
  my ($self, $out_file) = @_;
  my $max_length = 100000;
  
  my $file = path($out_file);
  
  my $data = $file->slurp;
  $data =~ s/\tnucleotide_match\t/\tmatch_part\t/gm;
  $data =~ s/\tprotein_match\t/\tmatch_part\t/gm;
  
  my @data = split(/\n/, $data);
  my %matches;
  
  my $new_data;
  
  foreach my $row (@data) {
    my ($region, $source, $start, $end, $strand, $name, $evalue) = $row =~
        /^(\S+)\t(\S+)\tmatch_part\t(\d+)\t(\d+)\t\S+\t(\S+).*Name=([^;]+).*evalue=([^;]+)/;
   
    if ($name) {
      my $id = join("_", $name, $region, $strand, $evalue);
      
      $matches{$id}{'region'} = $region;
      $matches{$id}{'source'} = $source;
      $matches{$id}{'strand'} = $strand;
      $matches{$id}{'name'}   = $name;
      if (!defined $matches{$id}{'start'} || $matches{$id}{'start'} > $start) {
        $matches{$id}{'start'} = $start;
      }
      if (!defined $matches{$id}{'end'} || $matches{$id}{'end'} < $end) {
        $matches{$id}{'end'} = $end;
      }
      $row .= ";Parent=$id";
      push @{$matches{$id}{'parts'}}, $row;
    } else {
      $new_data .= "$row\n";
    }
  }
  
  foreach my $id (keys %matches) {
    my $start = $matches{$id}{'start'};
    my $end   = $matches{$id}{'end'};
    my $name  = $matches{$id}{'name'};

    if ($end - $start < $max_length) {
      my $match = join("\t",
                       $matches{$id}{'region'},
                       $matches{$id}{'source'},
                       'match',
                       $start,
                       $end,
                       '.',
                       $matches{$id}{'strand'},
                       '.',
                       "ID=$id;Name=$name",
      );
      $new_data .= "$match\n";
      
      foreach my $part (@{$matches{$id}{'parts'}}) {
        $new_data .= "$part\n";
      }
    }
  }
  
  $file->spew($new_data);
}

sub Bio::EnsEMBL::Gene::summary_as_hash {
  my $self = shift;
  my %summary;
  
  $summary{'seq_region_name'} = $self->seq_region_name;
  $summary{'source'}          = $self->source;
  $summary{'start'}           = $self->seq_region_start;
  $summary{'end'}             = $self->seq_region_end;
  $summary{'strand'}          = $self->strand;
  $summary{'id'}              = $self->display_id;
  $summary{'Name'}            = $self->external_name;
  $summary{'biotype'}         = $self->biotype;
  $summary{'description'}     = $self->description;
  $summary{'version'}         = $self->version;
  
  return \%summary;
}

sub Bio::EnsEMBL::Transcript::summary_as_hash {
  my $self = shift;
  my %summary;
  
  my $parent_gene = $self->get_Gene();
  
  $summary{'seq_region_name'} = $self->seq_region_name;
  $summary{'source'}          = $parent_gene->source;
  $summary{'start'}           = $self->seq_region_start;
  $summary{'end'}             = $self->seq_region_end;
  $summary{'strand'}          = $self->strand;
  $summary{'id'}              = $self->display_id;
  $summary{'Parent'}          = $parent_gene->stable_id;
  $summary{'biotype'}         = $self->biotype;
  $summary{'description'}     = $self->description;
  $summary{'version'}         = $self->version;
  
  # Add xrefs
  my $xrefs = $self->get_all_xrefs();
  my (@db_xrefs, @go_xrefs);
  foreach my $xref (sort {$a->dbname cmp $b->dbname} @$xrefs) {
    my $dbname = $xref->dbname;
    if ($dbname eq 'GO') {
      push @go_xrefs, $xref->display_id;
    } else {
      $dbname =~ s/^RefSeq.*/RefSeq/;
      $dbname =~ s/^Uniprot.*/UniProtKB/;
      $dbname =~ s/^protein_id.*/NCBI_GP/;
      push @db_xrefs, "$dbname:".$xref->display_id;
    }
  }
  $summary{'Dbxref'} = \@db_xrefs if scalar(@db_xrefs);
  $summary{'Ontology_term'} = \@go_xrefs if scalar(@go_xrefs);
  
  # Check for seq-edits
  my $seq_edits = $self->get_all_SeqEdits();
  my @seq_edits;
  foreach my $seq_edit (@$seq_edits) {
    my ($start, $end, $alt_seq) =
      ($seq_edit->start, $seq_edit->end, $seq_edit->alt_seq);
    
    my $note = "This transcript's sequence has been ";
    if ($alt_seq eq '') {
      $note .= "annotated with a deletion between positions $start and $end";
    } elsif ($end < $start) {
      $note .= "annotated with an insertion before position $start: $alt_seq";
    } else {
      $note .= "replaced between positions $start and $end: $alt_seq";
    }
    
    push @seq_edits, $note;
  }
  $summary{'Note'} = \@seq_edits if scalar(@seq_edits);
  
  return \%summary;
}

sub Bio::EnsEMBL::Exon::summary_as_hash {
  my $self = shift;
  my %summary;
  
  $summary{'seq_region_name'} = $self->seq_region_name;
  $summary{'start'}           = $self->seq_region_start;
  $summary{'end'}             = $self->seq_region_end;
  $summary{'strand'}          = $self->strand;
  $summary{'id'}              = $self->display_id;
  $summary{'constitutive'}    = $self->is_constitutive;
  
  return \%summary;
}

sub Bio::EnsEMBL::RepeatFeature::summary_as_hash {
  my $self = shift;
  my %summary;
  
  $summary{'seq_region_name'} = $self->seq_region_name;
  $summary{'start'}           = $self->seq_region_start;
  $summary{'end'}             = $self->seq_region_end;
  $summary{'strand'}          = $self->strand;
  
  my $rc = $self->repeat_consensus;
  $summary{'Name'}  = $rc->name;
  $summary{'type'}  = $rc->repeat_type;
  $summary{'class'} = $rc->repeat_class;
  if ($rc->repeat_consensus =~ /^[^N]\S*/) {
    $summary{'repeat_consensus'} = $rc->repeat_consensus;
  }
  
  return \%summary;
}

sub Bio::EnsEMBL::BaseAlignFeature::summary_as_hash {
  my $self = shift;
  my %summary;
  
  $summary{'source'}              = $self->analysis->gff_source || $self->slice->source;
  $summary{'seq_region_name'}     = $self->seq_region_name;
  $summary{'start'}               = $self->seq_region_start;
  $summary{'end'}                 = $self->seq_region_end;
  $summary{'strand'}              = $self->strand;
  $summary{'id'}                  = undef;
  #$summary{'Parent'}              = $self->display_id.'_'.$self->seq_region_name.'_'.$self->strand.'_'.$self->p_value;
  $summary{'Name'}                = $self->display_id;
  $summary{'description'}         = $self->analysis->description;
  $summary{'score'}               = $self->score;
  $summary{'evalue'}              = $self->p_value;
  $summary{'percentage_identity'} = $self->percent_id;
  
  return \%summary;
}

1;
