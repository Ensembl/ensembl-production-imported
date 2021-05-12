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

package Bio::EnsEMBL::ENA::SRA::ExperimentAdaptor;

use strict;
use warnings;
use Carp;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use base qw( Bio::EnsEMBL::ENA::SRA::BaseSraAdaptor );
use Bio::EnsEMBL::ENA::SRA::Experiment;
use Data::Dumper;

sub get_by_accession {
    my ( $self, $acc ) = @_;
    return $self->get_by_accession_and_type($acc,'EXPERIMENT');
}


sub _hash_to_obj {
    my ( $self, $hash ) = @_;
    my $experiment = Bio::EnsEMBL::ENA::SRA::Experiment->new(
        -ACCESSION   => $hash->{accession},
        -CENTER_NAME => $hash->{center_name},
        -ALIAS       => $hash->{alias},
        -IDENTIFIERS => $hash->{IDENTIFIERS},
        -TITLE       => $hash->{TITLE},
        -DESIGN      => $hash->{DESIGN},
        -PLATFORM    => $hash->{PLATFORM},
        -LINKS       => $hash->{EXPERIMENT_LINKS}{EXPERIMENT_LINK},
        -ATTRIBUTES => $hash->{EXPERIMENT_ATTRIBUTES}{EXPERIMENT_ATTRIBUTE},
        -STUDY_ACCESSION => $hash->{STUDY_REF}{accession}
         );        
   return $experiment;
}
1;
