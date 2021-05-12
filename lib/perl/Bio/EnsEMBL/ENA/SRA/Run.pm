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

package Bio::EnsEMBL::ENA::SRA::Run;

use strict;
use warnings;
use Carp;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use base qw( Bio::EnsEMBL::ENA::SRA::BaseSraObject );
use Bio::EnsEMBL::ENA::SRA::BaseSraAdaptor qw(get_adaptor);
use Data::Dumper;

# 1 run belongs to 1 experiment (and hence to 1 study and multiple samples)
# 

sub new {
    my ( $proto, @args ) = @_;
    my $self = $proto->SUPER::new(@args);
    # accession, center_name, alias, identifiers hash
    (  $self->{title}, $self->{experiment_accession},
       )
      = rearrange( [ 'TITLE', 'EXPERIMENT_ACCESSION' ],
                   @args ) ;               
    return $self;
}

sub title {
    my ($self) = @_;
	return $self->{'title'};
}

sub experiment {
    my ($self) = @_;
    if(!defined $self->{experiment} && defined $self->{experiment_accession}) {    
        ($self->{experiment}) = @{get_adaptor('Experiment')->get_by_accession($self->{experiment_accession})};
    }
	return $self->{experiment};    
}
sub study { 
    my ($self) = @_;
    if(!defined $self->{study}) {
        $self->{studies} = $self->studies;
        warn "More than one study found for run " . $self->accession() if(scalar(@{$self->{studies}})>1);
        $self->{study} = $self->{studies}->[0];
    }
	return $self->{study};    
}
sub sample {
    my ($self) = @_;
    if(!defined $self->{sample}) {
        $self->{samples} = $self->samples;
        warn "More than one sample found for run ". $self->accession() if(scalar(@{$self->{samples}})>1);
        $self->{sample} = $self->{samples}->[0];
    }
	return $self->{sample};    
}
sub submission {
    my ($self) = @_;
    if(!defined $self->{submission}) {
        $self->{submissions} = $self->submissions;
        warn "More than one subission found for run". $self->accession() if(scalar(@{$self->{submissions}})>1);
        $self->{submission} = $self->{submissions}->[0];
    }
	return $self->{submission};    
}

1;
