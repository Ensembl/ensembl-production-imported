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

package Bio::EnsEMBL::EGPipeline::Xref::NamesAndDescriptions;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    'db_type' => 'core',
  };
}

sub write_output {
  my ($self) = @_;
  my $db_type = $self->param_required('db_type');
  my $timing  = $self->param_required('timing');

  my $dba = $self->get_DBAdaptor($db_type);
  my $ga  = $dba->get_adaptor('Gene');

  my %descriptions;
  my %gene_names;

  my $genes = $ga->fetch_all();
  foreach my $gene (@$genes) {
    if (defined $gene->description) {
      if ($gene->description =~ /Source:([^;]+)/) {
        $descriptions{$1}++;
      } else {
        $descriptions{'no source'}++;
      }
    } else {
      $descriptions{'null'}++;
    }

    if (defined $gene->external_db) {
      $gene_names{$gene->external_db}++;
    } else {
      $gene_names{'null'}++;
    }
  }
  
  $dba->dbc and $dba->dbc->disconnect_if_idle();

  $self->dataflow_output_id({}, 1);

  foreach my $db_name (keys %descriptions) {
    $self->dataflow_output_id(
    {
      'db_name' => $db_name,
      'total'   => $descriptions{$db_name},
      'timing'  => $timing,
    }, 2);
  }

  foreach my $db_name (keys %gene_names) {
    $self->dataflow_output_id(
    {
      'db_name' => $db_name,
      'total'   => $gene_names{$db_name},
      'timing'  => $timing,
    }, 3);
  }
}

1;
