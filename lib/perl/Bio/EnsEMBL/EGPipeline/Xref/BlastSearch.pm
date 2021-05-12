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

Bio::EnsEMBL::EGPipeline::Xref::BlastSearch

=head1 DESCRIPTION

Class for finding uniprot matches for a given sequence

=head1 Author

Dan Staines

=cut

package Bio::EnsEMBL::EGPipeline::Xref::BlastSearch;
use warnings;
use strict;
use Log::Log4perl qw/:easy/;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use LWP::UserAgent;
use XML::Simple;
use Data::Dumper;
use Carp;

=head1 CONSTRUCTOR
=head2 new


  Example    : $info = Bio::EnsEMBL::Utils::MetaData::GenomeInfo->new(...);
  Description: Creates a new info object
  Returntype : Bio::EnsEMBL::Utils::MetaData::GenomeInfo
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
  my ( $proto, @args ) = @_;
  my $class = ref($proto) || $proto;
  my $self = bless( {}, $class );
  $self->{logger} = get_logger();
  $self->{ua}     = LWP::UserAgent->new;
  return $self;
}

=head1 METHODS
=head2 logger
  Arg        : (optional) logger to set
  Description: Get logger
  Returntype : logger object reference
  Exceptions : none
  Caller     : internal
  Status     : Stable
=cut

sub logger {
  my ($self) = @_;
  if ( !defined $self->{logger} ) {
     $self->{logger} = get_logger();
  }
  return $self->{logger};
}

my $submit_url =
  'http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run/';
my $status_url =
  'http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/';
my $results_url =
  'http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/';

sub search {
  my ( $self, $seq, $collection, $stype, $program, $email ) = @_;

  $email //= $ENV{USER}.'@ebi.ac.uk',

  my $seq_hash = ref $seq eq 'HASH' ? $seq : { 1 => $seq };

  my %params = map {
    $_ => { email    => $email,
            sequence => $seq_hash->{$_},
            database => [$collection],
            program  => $program,
            exp      => '1e-4',
            stype    => $stype }
  } keys %$seq_hash;

  my $res = $self->run_blast( \%params );

  if (ref $seq eq 'HASH') {
     return $res;
  }
  return $res->{1};
} ## end sub search

sub run_blast {
  my ( $self, $inputs ) = @_;
  my $job_ids = {};
  while ( my ( $id, $params ) = each %$inputs ) {
     # submit job using post
     $job_ids->{$id} = $self->post( $submit_url, $params );
  }

  # poll until completed
  my $statuses = {};
  while () {
     while ( my ( $id, $job_id ) = each %$job_ids ) {
       if ( !defined $statuses->{$job_id} ) {
          my $status = $self->get( $status_url . $job_id );
          if ( $status ne 'RUNNING' ) {
            $statuses->{$job_id} = $status;
          }
       }
     }
     last if ( scalar( keys %$statuses ) == scalar( keys %$job_ids ) );
  }

  my $results = {};
  for my $id ( keys %{$inputs} ) {
     my $job_id = $job_ids->{$id};
     my $status = $statuses->{$job_id};
     if ( $status eq 'FINISHED' ) {
       # get results
       my $results_str = $self->get( $results_url . $job_id . '/xml' );
       # parse results as XML
       my $pres = XMLin($results_str)->{SequenceSimilaritySearchResult};
# service treats single and multiple hits differently in format, so need to check
       if ( $pres->{hits}{total} == 1 ) {
          $results->{$id} =
            { $pres->{hits}->{hit}->{id} => $pres->{hits}->{hit} };
       }
       else {
          $results->{$id} = $pres->{hits}->{hit};
       }

     }
     else {
       croak "BLAST completed with status $status";
     }
  }

  return $results;

} ## end sub run_blast

sub get {
  my ( $self, $url ) = @_;
  my $response = $self->{ua}->get($url);
  return $self->handle_response($response);
}

sub post {
  my ( $self, $url, $params ) = @_;
  my $response = $self->{ua}->post( $url, $params );
  return $self->handle_response($response);
}

sub handle_response {
  my ( $self, $response ) = @_;
  my $content = $response->content();
  if ( $response->is_error ) {
     my $error_message = '';
     # HTML response.
     if ( $content =~ m/<h1>([^<]+)<\/h1>/ ) {
       $error_message = $1;
     }
     #  XML response.
     elsif ( $content =~ m/<description>([^<]+)<\/description>/ ) {
       $error_message = $1;
     }
     croak 'Could not run BLAST: ' . $response->code .
       ' ' . $response->message . '  ' . $error_message;
  }
  return $content;

}

1;
