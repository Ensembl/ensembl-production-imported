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

package Bio::EnsEMBL::ENA::SRA::File;

use strict;
use warnings;
use Carp;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use base qw( Bio::EnsEMBL::Storable );
use Bio::EnsEMBL::ENA::SRA::BaseSraAdaptor qw(get_adaptor);
use LWP::Simple;
use HTTP::Status qw(:constants :is status_message);
use Data::Dumper;
use Bio::EnsEMBL::Utils::Exception qw/throw/;

sub new {
    my ( $proto, @args ) = @_;
    my $self = $proto->SUPER::new(@args);
    (  $self->{study_accession},
       $self->{sample_accession},
       $self->{experiment_accession},
       $self->{run_accession},
       $self->{file_name},
       $self->{file_size},
       $self->{url},
       $self->{md5},
       $self->{cols},

      ) = rearrange( [
          'STUDY_ACCESSION',      'SAMPLE_ACCESSION',
          'EXPERIMENT_ACCESSION', 'RUN_ACCESSION',
          'FILE_NAME',            'FILE_SIZE',
          'URL',                  'MD5', 'COLS'],
        @args );
    return $self;
}

sub cols {
    my ($self) = @_;
    return $self->{'cols'};
}

sub file_name {
    my ($self) = @_;
    return $self->{'file_name'};
}

sub file_size {
    my ($self) = @_;
    return $self->{'file_size'};
}

sub url {
    my ($self) = @_;
    return $self->{'url'};
}

sub md5 {
    my ($self) = @_;
    return $self->{'md5'};
}

sub experiment {
    my ($self) = @_;
    if ( !defined $self->{experiment} && defined $self->{experiment_accession} )
    {
        ( $self->{experiment} ) =
          @{ get_adaptor('Experiment')
              ->get_by_accession( $self->{experiment_accession} ) };
    }
    return $self->{experiment};
}

sub sample {
    my ($self) = @_;
    if ( !defined $self->{sample} && defined $self->{sample_accession} ) {
        ( $self->{sample} ) =
          @{ get_adaptor('Sample')
              ->get_by_accession( $self->{sample_accession} ) };
    }
    return $self->{sample};
}

sub study {
    my ($self) = @_;
    if ( !defined $self->{study} && defined $self->{study_accession} ) {
        ( $self->{study} ) =
          @{ get_adaptor('Study')->get_by_accession( $self->{study_accession} )
          };
    }
    return $self->{study};
}

sub run {
    my ($self) = @_;
    if ( !defined $self->{run} && defined $self->{run_accession} ) {
        ( $self->{run} ) =
          @{ get_adaptor('Run')->get_by_accession( $self->{run_accession} ) };
    }
    return $self->{run};
}
# Associations
#- associated with samples
#- associated with submissions
#- associated with ENA-FASTQ-FILES setstions

sub retrieve {
    my ($self, $dir) = @_;
    $dir ||= '.';
    my $file = $dir . '/' . $self->file_name();
    my $url = $self->url();
    
    if ($url =~ /^ftp\./) {

	# Won't work unless we add the domain protocol bit

	$url = "ftp://" . $url;
    }

    my $http_response_code = getstore($url, $file);

    if (! -f $file) {

	warn ("Response code, $http_response_code, " . status_message($http_response_code) . "\n");

	throw ("Failed to download file, $file from URL, $url!");
    }

    return $file;    
}

1;
