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

package Bio::EnsEMBL::ENA::SRA::BaseSraAdaptor;

use strict;
use warnings;
use Carp qw(cluck croak);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Data::Dumper;
use XML::Simple;
use LWP::Simple;
use Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyDBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw/throw/;
use Bio::EnsEMBL::ApiVersion;
use base qw(Exporter);

use constant PUBLIC_HOST => 'mysql-eg-publicsql.ebi.ac.uk';
use constant PUBLIC_PORT => 4157;
use constant PUBLIC_USER => 'anonymous';
use constant PUBLIC_PASS => '';

use Log::Log4perl qw(:easy);

my $logger = get_logger();

our @EXPORT = qw(get_adaptor taxonomy_adaptor get_linked_objects);

my $url = 'https://www.ebi.ac.uk/ena/browser/api/xml/%s';

my $adaptor_classes = {
    'experiment' => 'Bio::EnsEMBL::ENA::SRA::ExperimentAdaptor',
    'file'       => 'Bio::EnsEMBL::ENA::SRA::FileAdaptor',
    'run'        => 'Bio::EnsEMBL::ENA::SRA::RunAdaptor',
    'sample'     => 'Bio::EnsEMBL::ENA::SRA::SampleAdaptor',
    'study'      => 'Bio::EnsEMBL::ENA::SRA::StudyAdaptor',
    'submission' => 'Bio::EnsEMBL::ENA::SRA::SubmissionAdaptor',
};

my $adaptors = {};

sub get_adaptor {
    my ($type) = @_;
    $type = lc $type;
    my $adaptor = $adaptors->{$type};
    if ( !defined $adaptor ) {
        my $adaptor_class = $adaptor_classes->{$type};
        croak "No adaptor found for $type" unless defined $type;
        my $test_eval = eval "require $adaptor_class";
        if ( $@ or ( !$test_eval ) ) { croak($@) }
        $adaptor = $adaptor_class->new();
        $adaptors->{$type} = $adaptor;
    }
    return $adaptor;
}

sub new {
    my ( $proto, @args ) = @_;
    my $self = bless {}, $proto;
    ( $self->{taxonomy_adaptor}, $self->{ensembl_version} ) = rearrange( 'TAXONOMY_ADAPTOR', 'ENSEMBL_VERSION', @args );
    if ( !defined $self->{ensembl_version} ) {
      $self->{ensembl_version} = software_version();
    }
    if ( !defined $self->{taxonomy_adaptor} ) {
        my $dbname = 'ncbi_taxonomy_' . $self->{ensembl_version};
        my $tax_dba =
        Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyDBAdaptor->new(
          -host   => PUBLIC_HOST,
          -port   => PUBLIC_PORT,
          -user   => PUBLIC_USER,
          -pass   => PUBLIC_PASS,
          -dbname => $dbname,
        );
        $self->{taxonomy_adaptor} = $tax_dba->get_TaxonomyNodeAdaptor();
    }
    $self->clear_cache();
    return $self;
}

sub clear_cache {
    my ($self) = @_;
    $self->{cache} = {};
    return;
}

sub taxonomy_adaptor {
    my ($self) = @_;
    return $self->{taxonomy_adaptor};
}

sub _expand_acc {
    my ( $self, $acc_str ) = @_;
    my $accs = [];
    for my $acc ( @{ [$acc_str] } ) {
        for my $acc_elem ( split ',', $acc ) {
            if ( $acc_elem =~ m/[^-]+-[^-]+/ ) {
                my ( $start, $end ) = split '-', $acc_elem;
                $start =~ m/([A-Z]{3})([0-9]+)/;
                my $start_tag = $1;
                my $start_n   = 0 + $2;
                $end =~ m/([A-Z]{3})([0-9]+)/;
                my $end_tag = $1;
                croak "Could not parse $acc_elem" if ( $end_tag ne $start_tag );
                my $end_n = 0 + $2;

                for ( my $n = $start_n; $n <= $end_n; $n++ ) {
                    push @$accs, sprintf( "%s%06d", $start_tag, $n );
                }
            } else {
                push @$accs, $acc_elem;
            }
        }
    }
    return $accs;
} ## end sub _expand_acc

my $link_tags = {
    'Experiment' => 'ENA-EXPERIMENT',
    'File'       => 'ENA-FASTQ-FILES',
    'Run'        => 'ENA-RUN',
    'Sample'     => 'ENA-SAMPLE',
    'Study'      => 'ENA-STUDY',
    'Submission' => 'ENA-SUBMISSION',
};

sub get_linked_objects {
    my ( $obj, $link_type, $acc ) = @_;

    my $objs     = [];
    my $link_tag = $link_tags->{$link_type};
    croak "Could not find tag for link type $link_type"
      unless defined $link_tag;
    my $adaptor = get_adaptor($link_type);
    for my $link_hash ( @{ $obj->links() } ) {
        my $link = $link_hash->{XREF_LINK};
        if ( defined $link ) {
            if ( $link->{DB} eq $link_tag ) {
                for my $obj ( @{ $adaptor->get_by_accession( $link->{ID}, $acc ) } ) {
                    push @{$objs}, $obj;
                }
            }
        }
    }
    return $objs;
}

sub get_by_accession_and_type {
    my ( $self, $id, $key ) = @_;
    my $objs = [];
    for my $acc ( @{ $self->_expand_acc($id) } ) {
        my $obj = $self->{cache}{$acc};
        if ( !defined $obj ) {

	    $logger->debug("Getting object using url " . sprintf( $url, $acc ));

            my $doc = $self->get_document( sprintf( $url, $acc ) );
            if(!$doc->{$key}) {
            	throw("$acc does not appear to be of type $key");
            }
            $doc                 = $doc->{$key};
            $obj                 = $self->_hash_to_obj($doc);
            $obj->hash($doc);
            $self->{cache}{$acc} = $obj;
        }
        push @$objs, $obj;
    }

    return $objs;
}

sub get_document {
    my ( $self, $url ) = @_;
    my $content = get $url;
    croak "Couldn\'t get $url" unless $content;
    my $doc = XMLin($content);
    return $doc;
}

sub get_by_accession {
    croak "get_by_accession not implemented";
}

sub _hash_to_obj {
    croak "_hash_to_obj not implemented";
}

1;
