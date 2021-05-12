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

Bio::EnsEMBL::EGPipeline::RNAFeatures::EmailRNAReport

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::RNAFeatures::EmailRNAReport;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EmailReport');

use Email::MIME;
use Email::Sender::Simple;
use Path::Tiny;
use File::Path qw(make_path);

sub fetch_input {
  my ($self) = @_;
  my $run_cmscan         = $self->param_required('run_cmscan');
  my $run_trnascan       = $self->param_required('run_trnascan');
  my $load_mirbase       = $self->param_required('load_mirbase');
  my $cmscan_threshold   = $self->param_required('cmscan_threshold');
  my $trnascan_threshold = $self->param_required('trnascan_threshold');
  my $pipeline_dir       = $self->param_required('pipeline_dir');
  my $evalue_levels      = $self->param_required('evalue_levels');
  
  my $text = "The RNA Features pipeline has completed.\n";
  
  if ($run_cmscan) {
    my @evalue_levels = keys %$evalue_levels;
    foreach my $evalue_level (@evalue_levels) {
      if ($evalue_level > $cmscan_threshold) {
        delete $$evalue_levels{$evalue_level};
      }
    }
    @evalue_levels = keys %$evalue_levels;
    
    my $cmscan_file = "$pipeline_dir/cmscan.txt";
    $text .= "\nThe cmscan alignments are summarised in the table below, and in the attached plots.";
    $text .= $self->summarise($cmscan_file, 'E-value', \@evalue_levels);
    
  }
  
  if ($run_trnascan) {
    my $scores = [40];
    if ($trnascan_threshold < 40) {
      $scores = [0, 40];
    }
    
    my $trnascan_file = "$pipeline_dir/trnascan.txt";
    $text .= "\nThe trnascan alignments are summarised in the table below.";
    $text .= $self->summarise($trnascan_file, 'Score', $scores);
  }
  
  if ($load_mirbase) {
    my $mirbase_file = "$pipeline_dir/mirbase.txt";
    $text .= "\nThe mirbase load is summarised in the table below.";
    $text .= $self->summarise($mirbase_file, 'Score');
  }
  
  $self->param('text', $text);
}

sub run {
  my ($self) = @_;
  my $email         = $self->param_required('email');
  my $subject       = $self->param_required('subject');
  my $text          = $self->param_required('text');
  my $pipeline_dir  = $self->param_required('pipeline_dir');
  my $evalue_levels = $self->param_required('evalue_levels');

  # saving to $pipeline_dir
  if (defined $pipeline_dir && $text) {
    if (!-e $pipeline_dir) {
      make_path($pipeline_dir) or $self->throw("Failed to create directory '$pipeline_dir'");
    }
    my $path = "$pipeline_dir/rna_features_summary.txt";
    open( my $fh, ">", "$path")
      or die "Can't open > $path: $!";
    print $fh $text;
    close($fh);
  }
  
  # sending
  my @parts;
  push @parts, Email::MIME->create(body => $text);
  
  my $msg = Email::MIME->create(
    header => [
      From    => $email,
      To      => $email,
      Subject => $subject,
    ],
    parts => \@parts,
  );
  
  Email::Sender::Simple->send($msg);
}

sub summarise {
  my ($self, $file, $value_type, $thresholds) = @_;
  
  open(my $fh, '<', $file) or $self->throw("Failed to open $file: $!");
  
  my %summary;
  while (my $row = <$fh>) {
    my ($value, $species, $biotype, $name) = split(/\t/, $row);
    if (defined $thresholds) {
      foreach my $threshold (@$thresholds) {
        if ($value_type eq 'E-value') {
          if ($value <= $threshold) {
            $summary{$species}{$threshold}{$biotype}{$name}++;
          }
        } else {
          if ($value >= $threshold) {
            $summary{$species}{$threshold}{$biotype}{$name}++;
          }
        }
      }
    } else {
      $summary{$species}{''}{$biotype}{$name}++;
    }
  }
  
  close($fh);
  
  my @columns = ('Species', $value_type, '# Biotypes', '# Models', '# Alignments');
  my @results;
  foreach my $species (sort keys %summary) {
    foreach my $threshold (sort keys %{$summary{$species}}) {
      my ($biotypes, $models, $alignments) = (0, 0, 0);
      
      foreach my $biotype (keys %{$summary{$species}{$threshold}}) {
        $biotypes++;
        foreach my $name (keys %{$summary{$species}{$threshold}{$biotype}}) {
          $models++;
          $alignments += $summary{$species}{$threshold}{$biotype}{$name};
        }
      }
      
      push @results, [$species, $threshold, $biotypes, $models, $alignments];
    }
  }
  
  return $self->format_table('', \@columns, \@results);
}

1;
