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

Bio::EnsEMBL::Analysis::Runnable::DustMasker

=head1 SYNOPSIS

  my $dust = Bio::EnsEMBL::Analysis::Runnable::DustMasker->
  new(
      -query => $slice,
      -program => 'dustmasker',
     );
  $dust->run;
  my @repeats = @{$dust->output};

=head1 DESCRIPTION

Dust is a wrapper for the NCBI dustmasker program which runs the dust
algorithm to identify and mask simple repeats.

=cut

package Bio::EnsEMBL::Analysis::Runnable::DustMasker;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);

=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::DustMasker
  Arg [2]   : int, level
  Arg [3]   : int, linker
  Arg [4]   : int, window size
  Arg [5]   : int, split length
  Function  : create a new  Bio::EnsEMBL::Analysis::Runnable::DustMasker
  Returntype: Bio::EnsEMBL::Analysis::Runnable::DustMasker
  Exceptions: 
  Example   : 

=cut

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($level, $linker, $window_size, $split_length) =
    rearrange(['LEVEL', 'LINKER', 'WINDOW_SIZE', 'SPLIT_LENGTH',], @args);
  
  $self->program('dustmasker') if (!$self->program);
  $self->level($level) if ($level);
  $self->linker($linker) if ($linker);
  $self->window_size($window_size) if ($window_size);
  $split_length = 50000 unless defined $split_length;
  $self->split_length($split_length);
  
  return $self;
}

sub level {
  my $self = shift;
  $self->{'level'} = shift if (@_);
  return $self->{'level'};
}

sub linker {
  my $self = shift;
  $self->{'linker'} = shift if (@_);
  return $self->{'linker'};
}

sub window_size {
  my $self = shift;
  $self->{'window_size'} = shift if (@_);
  return $self->{'window_size'};
}

sub split_length {
  my $self = shift;
  $self->{'split_length'} = shift if (@_);
  return $self->{'split_length'};
}

=head2 run_analysis

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::DustMasker
  Arg [2]   : string, program name
  Function  : constructs a commandline and runs the program passed
  in, the generic method in Runnable isnt used as Dust doesnt
  fit this module
  Returntype: none
  Exceptions: throws if run failed because system doesnt
  return 0
  Example   : 

=cut

sub run_analysis {
  my ($self, $program) = @_;
  if (!$program) {
    $program = $self->program;
  }
  throw($program." is not executable Dust::run_analysis ") 
    unless($program && -x $program);
  my $command = $self->program;
  
  my $options = "";  
  $options .= " -level ".$self->level if ($self->level);
  $options .= " -linker ".$self->linker if ($self->linker);
  $options .= " -window ".$self->window_size if ($self->window_size);
  $self->options($options);
  
  $command .= "$options -in ".$self->queryfile." -out ".$self->resultsfile;
  
  system($command) == 0 or throw("FAILED to run ".$command);
}

=head2 parse_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::DustMasker
  Arg [2]   : string, filename
  Function  : open and parse the results file into repeat features
  Returntype: none 
  Exceptions: throws on failure to open or close output file
  Example   : 

=cut

sub parse_results{
  my ($self, $results) = @_;
  if(!$results){
    $results = $self->resultsfile;
  }
  my $ff = $self->feature_factory;
  my @output;
  open(DUST, $results) or throw("FAILED to open ".$results);
  LINE: while(<DUST>) {
    next LINE if (/^>/);
    if (/(\d+)\s\-\s(\d+)/) {
      my ($start, $end) = ($1, $2);
      # Dust results have 0-based start and end co-ordinates.
      $start++;
      $end++;
      
      my $rc = $ff->create_repeat_consensus('dust', 'dust', 'simple', 'N');
      my $rf = $ff->create_repeat_feature(
        $start, $end, 0, 0, $start, $end, $rc, $self->query->name, $self->query);
      
      if ($rf->length > $self->split_length) {
        my $converted_features = $self->convert_feature($rf);
        push(@output, @$converted_features) if ($converted_features);
      } else {
        push(@output, $rf);
      }
    }
  }
  $self->output(\@output);
  close(DUST) or throw("FAILED to close ".$results);
}

sub convert_feature {
  my ($self, $rf) = @_;
  
  my $ff = $self->feature_factory;
  my $projections = $rf->project('seqlevel');
  my @converted;
  my $feature_length = $rf->length;
  my $projected_length = 0;
  PROJECT:foreach my $projection (@$projections) {
    $projected_length += ($projection->from_end - $projection->from_start) +1;
  }
  
  my $percentage = 100;
  if($projected_length != 0) {
    $percentage = ($projected_length / $feature_length)*100;
  }
  if($percentage <= 75){
    return;
  }
  
  REPEAT:foreach my $projection (@$projections) {
    my $start = 1;
    my $end = $projection->to_Slice->length;
    my $slice = $projection->to_Slice;
    my $rc = $ff->create_repeat_consensus('dust', 'dust', 'simple', 'N');
    my $rf = $ff->create_repeat_feature($start, $end, 0, 0, $start,
                                        $end, $rc, $slice->name,
                                        $slice);
    my $transformed = $rf->transform($self->query->coord_system->name,
                                     $self->query->coord_system->version);
    if (!$transformed) {
      $self->throw("Failed to transform Dust region".
        $rf->seq_region_name.":".$rf->start."-".$rf->end);
    } else {
      $self->warning("Transformed Dust region".
        $rf->seq_region_name.":".$rf->start."-".$rf->end);
    }
    push(@converted, $transformed);
  }
  
  return \@converted;
}
