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

package Bio::EnsEMBL::ENA::SRA::Experiment;

use strict;
use warnings;
use Carp;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use base qw( Bio::EnsEMBL::ENA::SRA::BaseSraObject );
use Bio::EnsEMBL::ENA::SRA::BaseSraAdaptor qw(get_adaptor);
use Data::Dumper;

sub new {
    my ( $proto, @args ) = @_;
    my $self = $proto->SUPER::new(@args);
    # accession, center_name, alias, identifiers hash
    (  $self->{title}, $self->{design},
       $self->{platform}, $self->{study_accession} )
      = rearrange( [ 'TITLE',      
                     'DESIGN',
                     'PLATFORM', 'STUDY_ACCESSION'],
                   @args );
    return $self;
}

sub title {
    my ($self) = @_;
    return $self->{'title'};
}

sub design {
    my ($self) = @_;
    return $self->{'design'};
}

sub platform {
    my ($self) = @_;
    return $self->{'platform'};
}


sub study {
    my ($self) = @_;
    if(!defined $self->{study} && defined $self->{study_accession}) {    
        ($self->{study}) = @{get_adaptor('Study')->get_by_accession($self->{study_accession})};
    }
	return $self->{study};    
}

1;
