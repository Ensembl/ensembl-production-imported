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

package Bio::EnsEMBL::ENA::SRA::Submission;

use strict;
use warnings;
use Carp;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use base qw( Bio::EnsEMBL::ENA::SRA::BaseSraObject );
use Data::Dumper;

sub new {
    my ( $proto, @args ) = @_;
    my $self = $proto->SUPER::new(@args);
    # title, comment
    (  $self->{title}, $self->{comment}
       )
      = rearrange( [ 'TITLE', 'COMMENT' ],
                   @args ) ;               
    return $self;
}

sub title {
    my ($self) = @_;
	return $self->{'title'};
}
sub comment {
    my ($self) = @_;
	return $self->{'comment'};
}

1;
