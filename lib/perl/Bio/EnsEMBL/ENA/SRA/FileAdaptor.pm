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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.
 
=cut

package Bio::EnsEMBL::ENA::SRA::FileAdaptor;

use strict;
use warnings;
use Carp;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use base qw( Bio::EnsEMBL::ENA::SRA::BaseSraAdaptor );
use Bio::EnsEMBL::ENA::SRA::File;
use LWP::Simple;
use Data::Dumper;
use Bio::EnsEMBL::Utils::Exception qw/throw/;

use Log::Log4perl qw(:easy);

my $logger = get_logger();

sub get_by_accession {
  my ($self, $url, $acc) = @_;
  
  # As in the Experiment (or Run) xml the links are not correct, we correct here
  # Means it is now hardcoded, not ideal, but at least it is working
  my @fields = qw(
    run_accession
    fastq_ftp
    fastq_md5
    fastq_bytes
    library_layout
    secondary_study_accession
    secondary_sample_accession
    experiment_accession
    run_accession
  );
  my $fields_str = join(",", @fields);
  $url = "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$acc&result=read_run&fields=$fields_str";
  
  $logger->info("Get read files by accession, $acc, using url, '$url'...");

  return $self->parse_filelist($url);
}

sub parse_filelist {
  my ($self, $url) = @_;
  my $files = [];
  # get the document
  my $content = get $url;
  croak "Couldn\'t get content for, URL, $url" unless $content;
  
  # parse line by line
  open my $str_fh, '<', \$content;
  
  # First, get the columns names
  my $first_line = readline $str_fh;
  chomp $first_line;
  my @fields = split /\s+/, $first_line;
  
  while (my $line = <$str_fh>) {
    # Get each line
    chomp $line;
    my @cols = split '\t', $line;
    
    # Keep the column values in a hash
    my %vals;
    @vals{ @fields } = @cols;
    # Columns index .....
    if (
          ($vals{fastq_bytes} and $vals{fastq_bytes} =~ /;/)
          or
          ($vals{submitted_bytes} and $vals{submitted_bytes} =~ /;/)
    ) {
      if (defined $vals{library_layout} and $vals{library_layout} eq 'SINGLE') {
          $logger->warn('Library is annotated as SINGLE but it actually has several files: changed to PAIRED');
        }
        $vals{library_layout} = 'PAIRED';
    } else {
      if (defined $vals{library_layout} and $vals{library_layout} eq 'PAIRED') {
          $logger->warn('Library is annotated as PAIRED but it actually has one files: changed to SINGLE');
        }
        $vals{library_layout} = 'SINGLE';
    }
    
    if ($vals{library_layout} eq 'PAIRED') {

      # Both files are reported in the same column, so need to parse column to create file instances (separation is based on ';')
      $logger->info('PAIRED READS');

      # Two files then, so we need to be more clever
      my @filesizes;
      my @md5s;
      my @urls;
      if (defined $vals{fastq_bytes}) {
        @filesizes = split ';', $vals{fastq_bytes};
        @md5s      = split ';', $vals{fastq_md5};
        @urls      = split ';', $vals{fastq_ftp};
      }
      elsif (defined $vals{submitted_bytes}) {
        @filesizes = split ';', $vals{submitted_bytes};
        @md5s      = split ';', $vals{submitted_md5};
        @urls      = split ';', $vals{submitted_ftp};
      }
      $logger->info("Parsed ftp URLs for retrieving the fastq files, '" . join ("', '", @urls) . "'\n");

      # Parse the URLs to get the file names
      my @filenames = ();

      foreach my $url (@urls) {
        my @dirs = split '/', $url;
        my $filename = pop @dirs;
        push @filenames, $filename;
      } 

      my $i = 0;
      for my $url (@urls) {
        push @$files,
          Bio::EnsEMBL::ENA::SRA::File->new(
            -STUDY_ACCESSION      => $vals{secondary_study_accession},
            -SAMPLE_ACCESSION     => $vals{secondary_sample_accession},
            -EXPERIMENT_ACCESSION => $vals{experiment_accession},
            -RUN_ACCESSION        => $vals{run_accession},
            -FILE_NAME            => $filenames[$i],
            -FILE_SIZE            => $filesizes[$i],
            -MD5                  => $md5s[$i],
            -URL                  => $urls[$i],
            -COLS                 => [ map { $vals{ $_ } } @fields ],
          );
        $i++;
      }
    }
    
    else {
      $logger->info("SINGLE READS");
      
      # Get files metadata
      my $filesize = $vals{fastq_bytes};
      my $md5      = $vals{fastq_md5};
      my $url      = $vals{fastq_ftp};
      $filesize  ||= $vals{submitted_bytes};
      $md5       ||= $vals{submitted_md5};
      $url       ||= $vals{submitted_ftp};
      
      # can't find a filename, so use the URL as a filename
      # we could parse the url for getting the actual file name if necessary
      my @dirs = split '/', $url;
      my $filename = pop @dirs;
      $logger->info("Parsed ftp URL for retrieving the fastq file, '$url'\n");
      
      push @$files,
        Bio::EnsEMBL::ENA::SRA::File->new(
          -STUDY_ACCESSION      => $vals{secondary_study_accession},
          -SAMPLE_ACCESSION     => $vals{secondary_sample_accession},
          -EXPERIMENT_ACCESSION => $vals{experiment_accession},
          -RUN_ACCESSION        => $vals{run_accession},
          -FILE_NAME            => $filename,
          -FILE_SIZE            => $filesize,
          -MD5                  => $md5,
          -URL                  => $url,
          -COLS                 => [ map { $vals{ $_ } } @fields ],
        );
    }
  }
  return $files;
}

1;
