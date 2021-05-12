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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::tRNAscan

=head1 SYNOPSIS

  my $trnascan = Bio::EnsEMBL::Analysis::Runnable::tRNAscan->
  new(
      -query => $slice,
      -program => 'tRNAscan-SE',
     );
  $trnascan->run;
  my @features = @{$trnascan->output};

=head1 DESCRIPTION

tRNASCAN-SE predicts tRNA genes; this module runs it against a slice and
puts the results in the prediction_transcript table.

=cut

package Bio::EnsEMBL::Analysis::Runnable::tRNAscan;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);

=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::tRNAscan
  Arg [2]   : int, external_db_id: external_db_id from the external_db table
  Function  : create a new  Bio::EnsEMBL::Analysis::Runnable::tRNAscan
  Returntype: Bio::EnsEMBL::Analysis::Runnable::tRNAscan
  Exceptions: 
  Example   : 

=cut

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($pseudo, $threshold, $search_mode, $external_db_id) =
    rearrange(['PSEUDO', 'THRESHOLD', 'SEARCH_MODE', 'EXTERNAL_DB_ID',], @args);
  
  $self->program('tRNAscan-SE') if (!$self->program);
  $self->options(' -Q ') if (!$self->options);
  
  $pseudo = 0 unless $pseudo;
  $self->pseudo($pseudo);
  
  $self->threshold($threshold);
  
  $self->search_mode($search_mode);
  
  $self->external_db_id($external_db_id);
  
  return $self;
}

sub pseudo {
  my $self = shift;
  $self->{'pseudo'} = shift if (@_);
  return $self->{'pseudo'};
}

sub threshold {
  my $self = shift;
  $self->{'threshold'} = shift if (@_);
  return $self->{'threshold'};
}

sub search_mode {
  my $self = shift;
  $self->{'search_mode'} = shift if (@_);
  return $self->{'search_mode'};
}

sub external_db_id {
  my $self = shift;
  $self->{'external_db_id'} = shift if (@_);
  return $self->{'external_db_id'};
}

=head2 run_analysis

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::tRNAscan
  Arg [2]   : string, program name
  Function  : generates and executes the cmscan command
  Returntype: none
  Exceptions: throws if execution failed
  Example   : 

=cut

sub run_analysis {
  my ($self, $program) = @_;
  if (!$program) {
    $program = $self->program;
  }
  throw("$program is not executable") unless($program && -x $program);
  
  my $options = $self->options;
  unless ($self->pseudo) {
    $options .= " -D ";
  }
  if (defined $self->threshold) {
    $options .= " -X " . $self->threshold;
  }
  if ($self->search_mode) {
    if ($self->search_mode =~ /arch/i) {
      $options .= " -A ";
    } elsif ($self->search_mode =~ /bact/i) {
      $options .= " -B ";
    } elsif ($self->search_mode =~ /general/i) {
      $options .= " -G ";
    } elsif ($self->search_mode =~ /organ/i) {
      $options .= " -O ";
    }
  }
  $self->options($options);
  
  my $in_file  = $self->queryfile;
  my $out_file = $self->queryfile.'.out';
  $self->resultsfile($out_file);
  
  my $cmd = "$program $options -o $out_file $in_file";
  
  system($cmd) == 0 or throw("Failed to run: $cmd");
}

=head2 parse_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::CMScan
  Arg [2]   : string, filename
  Function  : open and parse the results file into sequence alignments
  Returntype: none 
  Exceptions: throws on failure to open or close output file
  Example   : 

=cut

sub parse_results {
  my ($self, $results) = @_;
  if (!$results) {
    $results = $self->resultsfile;
  }
  
  open(my $fh, $results) or throw("Failed to open results file '$results'");
  my $results_text = do { local $/; <$fh> };
  close($fh) or throw("Failed to close results file '$results'");
  
  $results_text =~ s/^Sequence.*\nName.*\n^\-+.*\n//m;
  my @lines = split(/\n/, $results_text);
  
  my @features;
  
  foreach my $line (@lines) {
    chomp $line;
    my (
      $seqname,
      undef,
      $start,
      $end,
      $aa_name,
      $anticodon,
      $intron_start,
      $intron_end,
      $score
    ) = split(/\s+/, $line);
  
    my $strand = 1;
    if ($start > $end) {
      ($start, $end) = ($end, $start);
      $strand = -1;
    }
    $anticodon =~ s/T/U/g;
    
    my ($length, $cigar);
    if ($intron_start == 0) {
      $length = $end - $start + 1;
      $cigar = $length.'M';
    } else {
      if ($intron_start > $intron_end) {
        ($intron_start, $intron_end) = ($intron_end, $intron_start);
      }
      my $exon_1_length = $intron_start - $start + 1;
      my $intron_length = $intron_end - $intron_start + 1;
      my $exon_2_length = $end - $intron_end + 1;
      $length = $exon_1_length + $exon_2_length;
      $cigar = $exon_1_length.'M'.$intron_length.'D'.$exon_2_length.'M';
    }
    
    my $rna_name = "tRNA-$aa_name";
    my $rna_desc = "$rna_name for anticodon $anticodon";

    my $biotype = "tRNA";
    $biotype .= "_pseudogene" if $aa_name eq 'Pseudo';
    
    my @attribs;
    push @attribs, $self->create_attrib('rna_gene_biotype', $biotype);
    push @attribs, $self->create_attrib('description',      $rna_desc);
    
    if (defined $score) {
      my $feature = Bio::EnsEMBL::DnaDnaAlignFeature->new(
        -slice          => $self->query,
        -start          => $start,
        -end            => $end,
        -strand         => $strand,
        -hstart         => 1,
        -hend           => $length,
        -hstrand        => 1,
        -score          => $score,
        -hseqname       => $rna_name,
        -cigar_string   => $cigar,
        -external_db_id => $self->external_db_id,
      );
      
      $feature->add_Attributes(@attribs);
      
      push (@features, $feature);
    }
  }
  
  $self->output(\@features);
}

sub create_attrib {
  my ($self, $key, $value) = @_;
  
  my $attrib = Bio::EnsEMBL::Attribute->new(
    -CODE  => $key,
    -VALUE => $value
  );
  
  return $attrib;
}

1;
