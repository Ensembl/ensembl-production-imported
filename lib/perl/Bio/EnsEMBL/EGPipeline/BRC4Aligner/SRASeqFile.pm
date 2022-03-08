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
use Digest::MD5;

sub param_defaults {
  my ($self) = @_;
  
  return {
    # Check the files downloaded (checksum)
    'check_fastq'  => 1,
  };
}

sub run {
  my ($self) = @_;
  
  my $work_dir = $self->param_required('work_dir');
  my $run_id   = $self->param_required('run_id');
  my $use_ncbi   = $self->param('use_ncbi');
  
  make_path($work_dir) unless -e $work_dir;
  
  my @runs; 
  if ($run_id =~ /^.RR/) {
    my $run_adaptor = get_adaptor('Run');
    @runs = @{$run_adaptor->get_by_accession($run_id)};
  } elsif ($run_id =~ /^.RS/) {
    my $sample_adaptor = get_adaptor('Sample');
    @runs = @{$sample_adaptor->get_by_accession($run_id)};
  } else {
    die "Unrecognized accession format '$run_id'";
  }
  die "No runs retrieved from accession '$run_id'" if @runs == 0;
  
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
  
  # Define paired with the number of files
  my $paired = 0;
  if (@files > 1) {
    $paired = 1;
  }
  elsif (@files == 1) {
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
        my $checksum_ok = $self->check_file($seq_file, $file->{md5});
        if (not $checksum_ok) {
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
        my $checksum_ok = $self->check_file($fastq_dl, $file->{md5});
        if (not $checksum_ok) {
          $self->throw("Downloaded file $fastq_dl is corrupted.\n");
        } else {
          rename $fastq_dl, $seq_file;
        }
        
        if (not -s $seq_file) {
          $self->throw("Retrieved file is empty '$fastq_dl'");
        }
      }
    } else {
      $self->throw("Cannot process file '$file_name'");
    }
  }

  return ($seq_file_1, $seq_file_2, $sam_file);
}

sub check_file {
  my ($self, $path, $md5) = @_;
  
  return 1 unless $self->param('check_fastq');
  
  print STDERR "Checking file $path...\n";
  
  open my $fileh, "<", $path or die("Couldn't read: $!");
  
  my $ctx = Digest::MD5->new;
  $ctx->addfile($fileh);
  
  my $dl_md5 = $ctx->hexdigest();
  
  print STDERR "Expected md5:   $md5\n";
  print STDERR "calculated md5: $dl_md5\n";
  
  return $dl_md5 eq $md5;
}

1;
