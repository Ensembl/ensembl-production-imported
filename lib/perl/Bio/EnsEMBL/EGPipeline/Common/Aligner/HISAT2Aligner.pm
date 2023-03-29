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

package Bio::EnsEMBL::EGPipeline::Common::Aligner::HISAT2Aligner;

use strict;
use warnings;
use base qw(Bio::EnsEMBL::EGPipeline::Common::Aligner);

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw);

use File::Basename;
use File::Spec::Functions qw(catdir);

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  (
    $self->{max_intron_length},
    $self->{gtf_file},
    $self->{strand_direction},
    $self->{is_paired},
    $self->{no_spliced},
    $self->{number_primary},
  ) =
  rearrange(
    [
      'MAX_INTRON_LENGTH', 'GTF_FILE', 'STRAND_DIRECTION', 'IS_PAIRED', 'NO_SPLICED', 'NUMBER_PRIMARY'
    ], @args
  );
  
  $self->{index_program} = 'hisat2-build';
  $self->{align_program} = 'hisat2';
  
  if ($self->{run_mode} eq 'lenient') {
    $self->{align_params} = ''; # Lenient params TBD
  } else {
    $self->{align_params} = '';
  }
  
  if ($self->{max_intron_length}) {
    $self->{align_params} .= " --max-intronlen $self->{max_intron_length} ";
  }

  if ($self->{strand_direction}) {
    my $strandness = "";
    if ($self->{strand_direction} eq 'forward') {
      if ($self->{is_paired}) {
        $strandness = "FR";
      } else {
        $strandness = "F";
      }
    } elsif ($self->{strand_direction} eq 'reverse') {
      if ($self->{is_paired}) {
        $strandness = "RF";
      } else {
        $strandness = "R";
      }
    } else {
      die "Unsupported strand_direction: " . $self->{strand_direction};
    }
    $self->{align_params} .= " --rna-strandness $strandness";
  }

  if ($self->{no_spliced}) {
    $self->{align_params} .= " --no-spliced-alignment";
  }

  if (defined($self->{number_primary})) {
    $self->{align_params} .= " -k " . $self->{number_primary};
  }

  if ($self->{aligner_dir}) {
    $self->{index_program} = catdir($self->{aligner_dir}, $self->{index_program});
    $self->{align_program} = catdir($self->{aligner_dir}, $self->{align_program});
  }
  
  return $self;
}

sub index_file {
  my ($self, $file) = @_;
  
  my $index_cmd  = $self->{index_program};
  my $index_name = $self->index_name($file);
  
  my $gtf_file = $self->{gtf_file};
  if ($gtf_file) {
    my $ss_file    = "$gtf_file.ss";
    my $exon_file  = "$gtf_file.exon";
    my $aligner_dir = $self->{aligner_dir};
    
    my $ss_cmd   = "python $aligner_dir/hisat2_extract_splice_sites.py $gtf_file > $ss_file";
    my $exon_cmd = "python $aligner_dir/hisat2_extract_exons.py $gtf_file > $exon_file";
    
    system($ss_cmd)   == 0 || throw "Cannot execute $ss_cmd: $@";
    system($exon_cmd) == 0 || throw "Cannot execute $exon_cmd: $@";
    $self->index_cmds($ss_cmd);
    $self->index_cmds($exon_cmd);
    
    $index_cmd .= " --ss $ss_file --exon $exon_file ";
  }
  
  $index_cmd .= " -p $self->{threads}";
  $index_cmd .= " $file $index_name";
  
  $self->run_cmd($index_cmd, 'index');
}

sub index_exists {
  my ($self, $file) = @_;
  
  my $index_name = $self->index_name($file);
  my $exists = -e "$index_name.1.ht2" ? 1 : 0;
  
  return $exists;
}

sub index_name {
  my ($self, $file) = @_;
  
  (my $index_name = $file) =~ s/\.\w+$//;
  if ($self->{gtf_file}) {
    $index_name .= '_tran';
  }
  return $index_name;
}

sub align {
  my ($self, $ref, $sam, $file1, $file2) = @_;
  
  my $index_name = $self->index_name($ref);
  if (defined $file2) {
   	$sam = $self->align_file($index_name, $sam, " -1 $file1 -2 $file2 ");
  } else {
   	$sam = $self->align_file($index_name, $sam, " -U $file1 ");
  }
  
  return $sam;
}

sub align_file {
  my ($self, $index_name, $sam, $files) = @_;
  
  my $format = $files =~ /fastq/ ? ' -q ' : ' -f ';
  
  my $align_cmd = $self->{align_program};
  $align_cmd   .= " --quiet ";
  $align_cmd   .= " -x $index_name ";
  $align_cmd   .= " -p $self->{threads} ";
  $align_cmd   .= " $self->{align_params} ";
  $align_cmd   .= " $format ";
  $align_cmd   .= " $files > $sam ";
  
  $self->run_cmd($align_cmd, 'align');
  
  return $sam;
}

1;
