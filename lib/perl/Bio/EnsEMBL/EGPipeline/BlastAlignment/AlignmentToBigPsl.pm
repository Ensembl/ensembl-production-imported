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

package Bio::EnsEMBL::EGPipeline::BlastAlignment::AlignmentToBigPsl;

use strict;
use warnings;
use JSON;
use Bio::EnsEMBL::Utils::Exception qw(throw); # warning

use File::Basename;
use File::Spec::Functions qw(catfile);

use base ('Bio::EnsEMBL::EGPipeline::FileDump::BaseDumper');

use constant _VALID_ALIGN_FEATURES => { map {$_=>1} qw/ProteinAlignFeature DnaAlignFeature/ };

use constant _BB_EXTRA_INDICES => qw/name url/;

use constant _URL_HIT_NAME => '<a href="/id/%s" rel="external">View hit in external database</a>';
# use constant _VB_URL_ALL_HITS => 'https://metazoa.ensembl.org/%s/Location/Genome?db=otherfeatures;ftype=%s;id=%s';
# https:/metazoa.ensembl.org/Aedes_aegypti_lvp/Location/Genome?db=otherfeatures;ftype=ProteinAlignFeature;id=_ALIGN_ID_


sub param_defaults {
  my ($self) = @_;
  return {
    %{$self->SUPER::param_defaults},
    db_type               => 'otherfeatures',
    feature_type          => ['ProteinAlignFeature'],
    join_align_feature    => 1,
    top_only_regions      => 1,
    default_only_coords   => 1,
    file_varname          => 'bigPsl_file',
    file_type             => 'bigPsl',
    auto_sql_file         => catfile(dirname(__FILE__), qw/.. Common Converters bigPsl.as/),
    swap_strands          => 0,
  };
}


sub run {
  my ($self) = @_;
  my $out_file            = $self->param('out_file');
  my $top_only_regions    = $self->param('top_only_regions');
  my $default_only_coords = $self->param('default_only_coords');

  my $species     = $self->param_required('species');
  my $db_type     = $self->param_required('db_type');
  my $logic_name  = $self->param_required('logic_name')->[0];

  my $feature_types      = $self->param_required('feature_type');
  my $join_align_feature = $self->param_required('join_align_feature');

  my $dba = $self->get_DBAdaptor($db_type);
  my $sla = $dba->get_adaptor('Slice');

  my %fta;
  foreach my $feature_type (@$feature_types) {
    if (exists _VALID_ALIGN_FEATURES->{$feature_type}) {
      $fta{$feature_type} = $dba->get_adaptor($feature_type);
    } else {
        $self->warning("unvalid feature type $feature_type (not ", join(", ", keys %{(_VALID_ALIGN_FEATURES)} ) ,")" );
    }
  }
  %fta or $self->throw("no valid feature types to use among @$feature_types");

  my $slices = $top_only_regions
    ? $sla->fetch_all('toplevel')
    : $sla->fetch_all();

#warn "toplevel slices number ", scalar(@$slices), "\n";

  $self->{written_features_count} = 0;
  open(my $out_fh, '>', $out_file) or $self->throw("Cannot open file $out_file: $!");

  foreach my $slice (@$slices) {
    next if ($default_only_coords && !$slice->coord_system()->is_default());

    foreach my $feature_type (@$feature_types) {
      my $features = $self->fetch_features($feature_type, $fta{$feature_type}, $logic_name, $slice);
      # warn "feature_type $feature_type features number ", scalar(@$features), "\n";
      $self->print_features($out_fh, $fta{$feature_type}, $features, $join_align_feature);
    }
  }

  close($out_fh);
  $self->warning(int($self->{written_features_count})." features in bigPsl file $out_file");
  $self->param('out_file', $out_file);
}

sub fetch_features {
  my ($self, $feature_type, $adaptor, $logic_name, $slice) = @_;
  
  my $features = $logic_name
    ? $adaptor->fetch_all_by_Slice($slice, $logic_name)
    : $adaptor->fetch_all_by_Slice($slice);

  return $features;
}

sub print_features {
  my ($self, $out_fh, $adaptor, $features, $join_align_feature) = @_;
  return if (!$features || !@{$features});

  # mb extract hit_ids and fetch them in separate loop
  my $features_as_dicts = [ grep {defined $_} map {$self->feature_to_dict($adaptor, $_)} grep {defined $_} @$features ];
  return unless (@$features_as_dicts);

  if ($join_align_feature) {
    $features_as_dicts = $self->join_features($features_as_dicts);
  }

  $features_as_dicts = [ sort by_pos_strand @$features_as_dicts ];
  foreach my $feature (@$features_as_dicts) {
    $self->print_feature($out_fh, $feature);
  }
}

sub by_pos_strand {
  $a->{chromStart} <=> $b->{chromStart}
                   ||
    $a->{chromEnd} <=> $b->{chromEnd}
                   ||
      $a->{strand} cmp $b->{strand}
}

sub by_hit_name_strand_pos {
        $a->{name} cmp $b->{name}
                   ||
      $a->{strand} cmp $b->{strand}
                   ||
  $a->{chromStart} <=> $b->{chromStart}
}


sub join_features {
  my ($self, $feats) = @_;
  my $out = [];

  return $out unless ($feats && @$feats);

  my $prev_id = '';
  my @parts = ();

  foreach my $feat (sort by_hit_name_strand_pos @$feats) {
    my $id = join("\t", $feat->{name}, $feat->{strand});
    if ($id ne $prev_id) {
      push @$out, @{ $self->combine_features(\@parts) } if @parts;
      $prev_id = $id;
      @parts = ();
    }
    push @parts, $feat;
  }
  push @$out, @{ $self->combine_features(\@parts) } if @parts;
  return $out;
}

sub combine_features {
  my ($self, $parts) = @_;
  return unless @$parts;

# do not merge by now
return $parts;

  my %stat = ();
  for my $it (@$parts) {
     my $pval = $it->{_p_value};
     $stat{$pval} = [] if (!exists $stat{$pval});
     push @{$stat{$pval}}, $it;
  }

  my @res = ();
  for my $pval (keys %stat) {
     my $pval_parts = $stat{$pval};
     my $out = { %{$pval_parts->[0]} };

     my $first_start = $out->{chromStart};
     my $o_first_start = $out->{oChromStart};

     $out->{blockCount}  = scalar(@$pval_parts);
     $out->{chromStarts} = join(',', map {$_->{chromStart} - $first_start} @$pval_parts);
     $out->{blockSizes}  = join(',', map {$_->{chromEnd} - $_->{chromStart}} @$pval_parts);
     $out->{chromEnd} = $pval_parts->[-1]->{chromEnd};

     $out->{oChromStarts} = join(',', map {$_->{oChromStart} - $o_first_start} @$pval_parts);
     $out->{oChromEnd} = $pval_parts->[-1]->{oChromEnd};
     #
     # score / etc / ???
     push @res, $out;
  }
  return \@res;
}

sub feature_to_dict {
  my ($self, $adaptor, $feat) = @_;

  if (! ( ref $feat && ( $feat->isa("Bio::EnsEMBL::DnaDnaAlignFeature")
                         || $feat->isa("Bio::EnsEMBL::DnaPepAlignFeature") ) ) ) {
       throw("feature must be a Bio::EnsEMBL::DnaDnaAlignFeature "
           . " or Bio::EnsEMBL::DnaPepAlignFeature, not a [". ref($feat). "]." );
  }

  my $slice = $feat->slice();
  if(!defined($slice) || !($slice->isa("Bio::EnsEMBL::Slice") or $slice->isa('Bio::EnsEMBL::LRGSlice')) ) {
    throw("A slice must be attached to the features to be stored.");
  }

  my $hseqname = $feat->hseqname();
  if ( !$hseqname ) {
    throw("DnaDnaAlignFeature must define an hseqname.");
  }

  my $hstart  = $feat->hstart();
  my $hend    = $feat->hend();
  my $hstrand = $feat->hstrand();

  $adaptor->_check_start_end_strand($hstart, $hend, $hstrand, $feat->hslice);

  my $cigar_string = $feat->cigar_string();
  if ( !$cigar_string ) {
    $cigar_string = $feat->length() . 'M';
    $self->warning( "*AlignFeature does not define a cigar_string.\n"
        . "Assuming ungapped block with cigar_line=$cigar_string ." );
  }

  # perc_id  = 100 * id / N
  #my $mismatches = int(100 * $feat->identical_matches * (1/$feat->percent_id -  1));
  #$mismatches = 0 if ($mismatches < 0);

  my $score = $feat->score;
  $score = 1000 if ($score > 1000);

  my $blocks = $self->cigar2blocks($cigar_string, $hstart, $hend);

  # swap strands to compact view for DnaPep Align alignments (xblast)
  my $do_strands_swap = $self->param('swap_strands'); 
  my $strand = $feat->strand;
  if ($feat->isa("Bio::EnsEMBL::DnaPepAlignFeature") && $strand < 0 && $do_strands_swap) {
    # swap strands, to compact view
    ($strand, $hstrand) = ($hstrand, $strand)
  }

  return {
    chrom       => $slice->seq_region_name,
    chromStart  => $feat->start-1,
    chromEnd    => $feat->end,

    name        => $feat->hseqname,
    score       => $score,
    strand      => $self->strand2char($strand),

    thickStart  => $feat->start-1,
    thickEnd    => $feat->start,
    reserved    => "255,215,0",

    blockCount  => join(",", $blocks->{'blockCount'}),
    blockSizes  => join(",", @{$blocks->{'blockSizes'}}),
    chromStarts => join(",", @{$blocks->{'chromStarts'}}), # relative to chromStart

    oChromStart  => $hstart,
    oChromEnd    => $hend,
    oStrand      => $self->strand2char($hstrand),
    oChromSize   => $hend, 
    oChromStarts => join(",", @{$blocks->{'oChromStarts'}}), # relative to oChromStart

    oSequence    => $cigar_string,
    oCDS         => $self->browser_url(_URL_HIT_NAME, $feat->hseqname) . "pVal: ". $feat->p_value,

    chromSize    => $slice->seq_region_length(), 

    match        => ($feat->identical_matches() or 0),
    misMatch     => 0,
    repMatch     => 0,
    nCount       => 0,
    
    seqType      => 2, # amino acid

    # used for merging
    _p_value    => $feat->p_value,
    _percent_id => $feat->percent_id,
    _hcoverage  => $feat->hcoverage,
  };
}

sub cigar2blocks {
  my ($self, $cigar_raw, $h_start, $h_end) = @_;

  my $cigar = $cigar_raw;
  $cigar =~ s/(\d+)(\D+)/$2 $1,/g;
  $cigar =~ s/,$//;
  
  my @chromStarts = ();
  my @blockSizes  = ();
  my @oChromStarts = ();

  my $start = 0;
  my $hstart = 0;
  for my $pair (split /,/, $cigar) {
    my ($c, $count) = split / /, $pair;
    if ($c eq "M") {
      push @chromStarts, $start;  
      push @oChromStarts, $hstart;  
      push @blockSizes, $count;  
      $start += $count;
      $hstart += int($count / 3);  
      next;
    }
    if ($c eq "I") {
      $start += $count;
      next;
    }
    if ($c eq "D") {
      push @chromStarts, $start;  
      push @oChromStarts, $hstart;  
      $hstart += int($count / 3);  
      push @blockSizes, 0;  
      next;
    }
  } 

  return {
       blockCount    => scalar(@chromStarts), 
       chromStarts   => \@chromStarts,
       blockSizes    => \@blockSizes,
       oChromStarts  => \@oChromStarts,
     };
}

sub print_feature {
  my ($self, $out_fh, $feat) = @_;

  my @feat_cols = map {!defined $feat->{$_} && '.' || $feat->{$_}} @{ $self->get_fields_order() };
  ++$self->{written_features_count};
  print $out_fh join("\t", map {!defined $_ && '.' || $_} @feat_cols), "\n";
}

sub update_fields_order {
  my ($self) = @_;

  if (defined $self->{_fields_in_order}) {
    return $self->{_fields_in_order};
  }

  my @TEMPLATE = ();
  my $auto_sql_file = $self->param('auto_sql_file');
  open(my $as_fh, "<", $auto_sql_file)
    or die("Can't open auto_sql_file:  $auto_sql_file. $!");
  while (<$as_fh>) {
    chomp;
    next if m/^\s*$/;
    push @TEMPLATE, $_;
  }
  close($as_fh);

  my $fields_in_order = [
      map { m/^\s*+\w+(?:\[\s*\w+\s*\])?\s+(\w+)\s*;/ && $1 }
        grep /;/, @TEMPLATE
  ];
  
  warn "ORDER OF bigPsl FIELDS: ", join(", ", @$fields_in_order);

  $self->{_fields_in_order} = $fields_in_order;
}

sub get_fields_order {
  my ($self) = @_;
  if (!defined $self->{_fields_in_order}) {
    $self->update_fields_order();
  }
  return $self->{_fields_in_order};
}

sub write_output {
  my ($self) = @_;

  my $dataflow_output = {
    out_file => $self->param('out_file'),
  };
  
  $self->dataflow_output_id($dataflow_output, 1);
}

sub browser_url {
  my ($self, $template, @params) = @_;
  #return ".";
  # validation
  my $url = sprintf($template, @params);
  return $url;
}

sub strand2char {
  my ($self, $strand) = @_;
  return $strand == 1
    ? '+'
    : $strand == -1
      ? '-'
      : '.'
  ;
}

1;
