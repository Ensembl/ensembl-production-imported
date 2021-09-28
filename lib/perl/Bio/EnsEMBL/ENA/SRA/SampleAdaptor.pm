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

package Bio::EnsEMBL::ENA::SRA::SampleAdaptor;

use strict;
use warnings;
use Carp;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use base qw( Bio::EnsEMBL::ENA::SRA::BaseSraAdaptor );
use Bio::EnsEMBL::ENA::SRA::Sample;
use Data::Dumper;

my $taxurl = 
'https://www.ebi.ac.uk/ena/browser/api/xml/Taxon:%s';

sub new {
    my ( $proto, @args ) = @_;
    my $self = $proto->SUPER::new(@args);
}

sub get_by_accession {
    my ( $self, $acc ) = @_;
    return $self->get_by_accession_and_type($acc,'SAMPLE');
}

sub _hash_to_obj {
    my ( $self, $hash ) = @_;
    my $taxid = $hash->{SAMPLE_NAME}{TAXON_ID};
    if(!defined $taxid) {
        confess Dumper($hash);
    }
    my $taxon = $self->taxonomy_adaptor()
                       ->fetch_by_taxon_id( $taxid);
    my $sample =
      Bio::EnsEMBL::ENA::SRA::Sample->new(
                     -ACCESSION   => $hash->{accession},
                     -CENTER_NAME => $hash->{center_name},
                     -ALIAS       => $hash->{alias},
                     -IDENTIFIERS => $hash->{IDENTIFIERS},
                     -TITLE       => $hash->{TITLE},
                     -DESCRIPTION       => $hash->{DESCRIPTION},
                     -TAXON       => $taxon,
                     -LINKS      => $hash->{SAMPLE_LINKS}{SAMPLE_LINK},
                     -ATTRIBUTES => $hash->{SAMPLE_ATTRIBUTES}{SAMPLE_ATTRIBUTE}
      );

    return $sample;
} ## end sub _hash_to_sample

sub get_by_taxon_id {
    my ( $self, $taxid ) = @_;
    my $samples = [];
    for my $accession ( @{ $self->get_accessions_for_taxon_id($taxid) } ) {
        $samples = [@{$self->get_by_accession($accession)},@$samples];
    }
    return $samples;
}

sub get_accessions_for_taxon_id {
    my ( $self, $taxid ) = @_;
    my $samples = $self->get_document(sprintf($taxurl,$taxid)
    );
    if ( !defined $samples->{SAMPLE} ) {
        $samples->{SAMPLE} = [];
    } elsif ( ref( $samples->{SAMPLE} ) eq 'HASH' ) {
        $samples->{SAMPLE} = [ $samples->{SAMPLE} ];
    }
    my $accs = [];
    for my $sample (@{$samples->{SAMPLE}}) {
        push @$accs,$sample->{accession};
    }
    return $accs;
}
1;
