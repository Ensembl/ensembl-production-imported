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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::SRASeqFileFromNCBI;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::EnsEMBL::ENA::SRA::BaseSraAdaptor qw(get_adaptor);

use File::Path qw(make_path);
use File::Spec::Functions qw(catdir);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use LWP::Simple qw(getstore);
use Carp;
use Cwd qw(getcwd);

sub run {
  my ($self) = @_;
  
  my $work_dir = $self->param_required('work_dir');
  my $run_id   = $self->param_required('run_id');
  
  make_path($work_dir) unless -e $work_dir;
  
  my $run_adaptor = get_adaptor('Run');
  my @runs = @{$run_adaptor->get_by_accession($run_id)};
  
  my @output = ();
  
  foreach my $run (@runs) {
    my ($seq_file_1, $seq_file_2, $sam_file) = $self->retrieve_files($work_dir, $run);
    
    push @output,
      {
        'run_seq_file_1' => $seq_file_1,
        'run_seq_file_2' => $seq_file_2,
        'sam_file'   => $sam_file,
      };
  }
  
  $self->param('output', \@output);
}

sub write_output {
  my ($self) = @_;
  
  foreach my $output (@{$self->param('output')}) {
    $self->dataflow_output_id($output, 2);
  }
}

sub retrieve_files {
  my ($self, $work_dir, $run) = @_;
  
  my $sra_dir = $self->param('sra_dir');
  
  my $run_acc = $run->accession;
  my $sam_file = catdir($work_dir, "$run_acc.sam");
  
  # Try to use existing file (previously downloaded or manually added)
  {
      my ($seq_file_1, $seq_file_2);
      # Try to use a single file
      $seq_file_1 = catdir($work_dir, "$run_acc\_all.fastq.gz");
      if (-s $seq_file_1) {
          return ($seq_file_1, $seq_file_2, $sam_file);
      }
      # Or try to use a pair
      $seq_file_1 = catdir($work_dir, "$run_acc\_all_1.fastq.gz");
      $seq_file_2 = catdir($work_dir, "$run_acc\_all_2.fastq.gz");
      if (-s $seq_file_1 and -s $seq_file_2) {
          return ($seq_file_1, $seq_file_2, $sam_file);
      }
  }
  my ($seq_file_1, $seq_file_2);
  
  my $sra_acc = $run_acc;
  my $download_sra = 0;
  
  my $previous_dir = getcwd();
  chdir $work_dir;
  if ($download_sra) {
    # Retrieve the SRA file
    my $SRA_ROOT = "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/%s/%s/%s/%s.sra";
    my @sra_pars = (
      substr($run_acc, 0, 3),
      substr($run_acc, 0, 6),
      $run_acc,
      $run_acc
    );
    my $sra_url = sprintf($SRA_ROOT, @sra_pars);

    # Download it
    # We HAVE to be in the directory where the files will be extracted...
    my $sra_file = "$run_acc.sra";
    my $http_response_code = getstore($sra_url, $sra_file);
    die "Could not download file $sra_file from $sra_url" if not -f $sra_file;
    $sra_acc = $sra_file;
  }
  
  # Next, extract the fastq file(s)
  warn "Getting file from fastq-dump";
  
  my $sra_dump_bin = $sra_dir ? "$sra_dir/fastq-dump": "fastq-dump";
  
  # NB: -B ensures we get base-space data, not color-space
  my $sra_command = "$sra_dump_bin -I -B --split-files $sra_acc";
  print("SRA download command: '$sra_command'\n");
  my $stats = `$sra_command`;
  
  # The downloaded files should be named {$run_acc}_1.gz or {$run_acc}.gz
  # Check file names and rename to our standard names
  
  my $lone_fastq = "$run_acc.fastq";
  my $downloaded_1 = "${run_acc}_1.fastq";
  my $downloaded_2 = "${run_acc}_2.fastq";

  if (-f $lone_fastq) {
    $seq_file_1 = "${run_acc}_all.fastq";
    rename $lone_fastq, $seq_file_1;

    `gzip $seq_file_1`;
    $seq_file_1 .= '.gz';
  } elsif (-f $downloaded_1 and -f $downloaded_2) {
    # Check both files have the same number of lines
    my $lines1 = `cat $downloaded_1 | wc -l`; chomp $lines1;
    my $lines2 = `cat $downloaded_2 | wc -l`; chomp $lines2;
    if ($lines1 != $lines2) {
      die("Files have different line counts ($lines1 vs $lines2): $downloaded_1 and $downloaded_2");
    }

    $seq_file_1 = "${run_acc}_all_1.fastq";
    $seq_file_2 = "${run_acc}_all_2.fastq";
    rename $downloaded_1, $seq_file_1;
    rename $downloaded_2, $seq_file_2;
    `gzip $seq_file_1`;
    $seq_file_1 .= '.gz';
    `gzip $seq_file_2`;
    $seq_file_2 .= '.gz';

  } elsif (-f $downloaded_1) {
    warn("Using $downloaded_1 as the sole file for $run_acc");
    $seq_file_1 = "${run_acc}_all.fastq";
    rename $downloaded_1, $seq_file_1;
    `gzip $seq_file_1`;
    $seq_file_1 .= '.gz';
  } else {
    chdir($previous_dir);
    die "Can't find the fastq files extracted for $run_acc downloaded in $work_dir";
  }
  
  # Remove sra file
  unlink $sra_acc if -s $sra_acc;

  # Back to the root dir: need to update the seq files path
  chdir($previous_dir);
  $seq_file_1 = catdir($work_dir, $seq_file_1);
  $seq_file_2 = catdir($work_dir, $seq_file_2) if $seq_file_2;
  
  return ($seq_file_1, $seq_file_2, $sam_file);
}

1;
