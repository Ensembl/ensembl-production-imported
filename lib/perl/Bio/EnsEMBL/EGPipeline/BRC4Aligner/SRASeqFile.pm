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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::SRASeqFile;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::EnsEMBL::ENA::SRA::BaseSraAdaptor qw(get_adaptor);

use File::Path qw(make_path);
use File::Spec::Functions qw(catdir);
use IO::Zlib;
use Try::Tiny;

sub param_defaults {
  my ($self) = @_;
  
  return {
    # Check the files downloaded: thorough but takes a long time
    'check_fastq'  => 0,
  };
}

sub run {
  my ($self) = @_;
  
  my $work_dir = $self->param_required('work_dir');
  my $run_id   = $self->param_required('run_id');
  my $use_ncbi   = $self->param('use_ncbi');
  
  make_path($work_dir) unless -e $work_dir;
  
  my $run_adaptor = get_adaptor('Run');
  my @runs = @{$run_adaptor->get_by_accession($run_id)};
  
  my @output = ();
  my @output_failed = ();
  
  foreach my $run (@runs) {
    try {
      my ($seq_file_1, $seq_file_2, $sam_file) = $self->retrieve_files($work_dir, $run);

      push @output,
      {
        'run_seq_file_1' => $seq_file_1,
        'run_seq_file_2' => $seq_file_2,
      };
    } catch {
      my $error = $_;

      if ($use_ncbi) {
        print STDERR "ERROR: $error\n";
        push @output_failed, { run_id => $run_id };
      } else {
        my $run_id = $run->accession;
        die "Could not download the files for $run_id: $error\n";
      }
    };
  }
  
  $self->param('output', \@output);
  $self->param('output_failed', \@output_failed);
}

sub write_output {
  my ($self) = @_;
  
  foreach my $output (@{$self->param('output')}) {
    $self->dataflow_output_id($output, 2);
  }

  foreach my $output (@{$self->param('output_failed')}) {
    $self->dataflow_output_id($output, 3);
  }
}

sub retrieve_files {
  my ($self, $work_dir, $run) = @_;
  
  my $run_acc = $run->accession;
  my @files = grep { $_->file_name() } @{$run->files()};
  
  my $experiment = $run->experiment();
  my $paired = defined $experiment->design()->{LIBRARY_DESCRIPTOR}{LIBRARY_LAYOUT}{PAIRED};
  
  # Single but several files?? Treat as single, as it is likely an error in SRA
  if (not $paired and @files > 1) {
    warn "Experiment is SINGLE, but there are several files. Changed to PAIRED.\n";
    $paired = 1;
  }
  if ($paired and @files == 1) {
    warn "Experiment is PAIRED, but there is only one file. Changed to SINGLE.\n";
    $paired = 0;
  } elsif (@files == 0) {
    die("There is no file to download for $run_acc");
  }
  print STDERR "There are ".scalar(@files)." files to download from $run_acc\n";
  print STDERR join(", ", map { $_->file_name() } @files)."\n";
  
  my ($seq_file_1, $seq_file_2) = ("", "");
  my $sam_file = catdir($work_dir, "$run_acc.sam");
  
  # Prepare files names
  if ($paired) {
    $seq_file_1 = catdir($work_dir, "$run_acc\_all_1.fastq.gz");
    $seq_file_2 = catdir($work_dir, "$run_acc\_all_2.fastq.gz");
  } else {
    $seq_file_1 = catdir($work_dir, "$run_acc\_all.fastq.gz");
  }
  
  my $new_files = 0;
  my %counts = ();

  # Retrieve each file
  FILE: for my $file (@files) {
    my $file_name = $file->file_name();
    print STDERR "Process file $file_name\n";
    
    # Decide the file name to use
    if ($file_name =~ /\.fastq/) {
      my $seq_file = $seq_file_1;
      
      # Choose a name for a pair
      if ($paired) {
        if ($file_name =~ /_1.fastq/) {
          $seq_file = $seq_file_1;
        } elsif ($file_name =~ /_2.fastq/) {
          $seq_file = $seq_file_2;
        } else {
          warn("Expected paired-end files, but got unnumbered file: '$file_name'\n");
          next FILE;
        }
      }
      
      # Reuse files if possible
      my $fastq_final = $seq_file;
      my $fastq_dl = catdir($work_dir, $file_name);

      # Reuse fastq file if it was succesfully unzipped
      my $reuse_file = 0;
      if (-s $seq_file and not -s $fastq_dl) {
        # Check the file first
        my $seq_count = $self->check_file($seq_file);
        $counts{$seq_file} = $seq_count;
        if (not $seq_count) {
          print STDERR "Current file $seq_file is corrupted.\n";
          $reuse_file = 0;
        } else {
          $reuse_file = 1;
        }
      }
      
      if ($reuse_file) {
        print STDERR "File $seq_file is already there. Skipping.\n";
        next FILE;
      } else {
        $new_files = 1;
        print STDERR "Downloading $fastq_dl for $seq_file.\n";
        unlink ($seq_file, $fastq_dl);
        $file->retrieve($work_dir);

        # Check downloaded file
        my $seq_count = $self->check_file($fastq_dl);
        if (not $seq_count) {
          $self->throw("Downloaded file $fastq_dl is corrupted.\n");
        } else {
          rename $fastq_dl, $seq_file;
          $counts{$seq_file} = $seq_count;
        }
        
        if (not -s $seq_file) {
          $self->throw("Retrieved file is empty '$fastq_dl'");
        }
      }
    } else {
      $self->throw("Cannot process file '$file_name'");
    }
  }

  # Check for files corruption
  if ($new_files) {
    if ($paired) {
      my $reads1 = $counts{$seq_file_1};
      my $reads2 = $counts{$seq_file_2};
      
      die("File $seq_file_1 is likely corrupted") if not $reads1;
      die("File $seq_file_2 is likely corrupted") if not $reads2;
      die("Files $seq_file_1 and $seq_file_2 have different read count: $reads1 vs $reads2") if $reads1 != $reads2;
    } else {
      my $reads = $counts{$seq_file_1};
      die("File $seq_file_1 is likely corrupted") if not $reads;
    }
  }
  
  return ($seq_file_1, $seq_file_2, $sam_file);
}

sub check_file {
  my ($self, $path) = @_;
  
  return 1 unless $self->param('check_fastq');
  
  print STDERR "Checking file $path...\n";
  
  # Check that all reads are whole, and return a count
  # If it is not whole, return undef
  
  my $read_count = 0;
  my $whole_count = 0;
  
  my $fh = IO::Zlib->new($path, "rb");
  if ($fh) {
    while (my $read = readline $fh) {
      # Read lines 4 by 4
      my $seq = readline $fh;
      my $strand = readline $fh;
      my $quality = readline $fh;
      if ($read =~ /^@/ and $seq and $quality and length($seq) == length($quality)) {
        $whole_count++;
      }
      $read_count++;
    }
    $fh->close;
  }
  
  if ($read_count == $whole_count) {
    return $read_count;
  } else {
    return;
  }
}

1;
