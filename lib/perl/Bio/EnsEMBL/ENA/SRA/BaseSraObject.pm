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

package Bio::EnsEMBL::ENA::SRA::BaseSraObject;

use strict;
use warnings;
use Carp;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use base qw( Bio::EnsEMBL::Storable );
use Bio::EnsEMBL::ENA::SRA::BaseSraAdaptor qw(get_linked_objects);

sub new {
    my ( $proto, @args ) = @_;
    my $self = $proto->SUPER::new(@args);
    # accession, center_name, alias, identifiers hash
    (  $self->{accession},   $self->{center_name}, $self->{alias},
       $self->{identifiers}, $self->{links},       $self->{attributes},$self->{hash}  )
      = rearrange( [ 'ACCESSION', 'CENTER_NAME', 'ALIAS', 'IDENTIFIERS',
                     'LINKS',     'ATTRIBUTES', 'HASH' ],
                   @args );
    $self->{links}       ||= [];
    $self->{attributes}  ||= [];
    $self->{identifiers} ||= [];
    return $self;
}

sub links {
    my ($self) = @_;
    return $self->{'links'};
}

sub attributes {
    my ($self) = @_;
    return $self->{'attributes'};
}

sub accession {
    my ($self) = @_;
    return $self->{'accession'};
}

sub center_name {
    my ($self) = @_;
    return $self->{'center_name'};
}

sub alias {
    my ($self) = @_;
    return $self->{'alias'};
}

sub identifiers {
    my ($self) = @_;
    return $self->{'identifiers'};
}

sub hash {
    my ($self,$hash) = @_;
    if(defined $hash) {
        $self->{hash} = $hash;
    }
    return $self->{'hash'};
}



## Functions to follow links to various types of SRA object. Each
## function returns a list of objects of the requested type. If no
## links exist between the current object type and the requested
## object type, an empty list is returned.

## TODO: Currently no attempt is made to recurse through object links
## to find objects of the requested type.

## TODO: No attempt is made to sanity check the various links.

sub experiments {
    my ($self) = @_;
    if(!defined $self->{experiments}) {
        $self->{experiments} = get_linked_objects($self, 'Experiment');
    }
    return $self->{experiments};
}

sub files {
     my ($self) = @_;
     if(!defined $self->{files}) {
        $self->{files}       = get_linked_objects($self, 'File', $self->accession());
     }
     return $self->{files};
}

sub runs {
    my ($self) = @_;
    if(!defined $self->{runs}) {
        $self->{runs}        = get_linked_objects($self, 'Run');
    }
    return $self->{runs};
}

sub samples {
    my ($self) = @_;
    if(!defined $self->{samples}) {
        $self->{samples}     = get_linked_objects($self, 'Sample');
    }
    return $self->{samples};
}

sub studies {
    my ($self) = @_;
    if(!defined $self->{studies}) {
        $self->{studies}     = get_linked_objects($self, 'Study');
    }
    return $self->{studies};
}

sub submissions {
    my ($self) = @_;
    if(!defined $self->{submissions}) {
        $self->{submissions} = get_linked_objects($self, 'Submission');
    }
	return $self->{submissions};
}

1;
