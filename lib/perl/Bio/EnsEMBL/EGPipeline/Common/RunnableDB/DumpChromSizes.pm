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

package Bio::EnsEMBL::EGPipeline::Common::RunnableDB::DumpChromSizes;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::FileDump::BaseDumper');
use File::Spec::Functions qw/splitpath/;

sub param_defaults {
  my ($self) = @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'top_only_regions'      => 1,
    'default_only_coords'   => 1,
    'file_type'             => 'tsv',
  };
}

sub run {
  my ($self) = @_;
  
  my $out_file            = $self->param('out_file');
  my $top_only_regions    = $self->param('top_only_regions');
  my $default_only_coords = $self->param('default_only_coords');

  my $dba = $self->core_dba();

  my %cs_interest_id ={};
  my $csa = $dba->get_adaptor('CoordSystem');
  my $sla = $dba->get_adaptor('Slice');

  my $slices = $top_only_regions
    ? $sla->fetch_all('toplevel')
    : $sla->fetch_all();

  my ( $_vname, $output_dir, $_fname) = splitpath($out_file);
  mkdir $output_dir unless -e $output_dir;
  open(my $fh, '>', $out_file) or $self->throw("Cannot open file $out_file: $!");

  for my $sl ( @{ $slices } ) {
    next if ($default_only_coords && !$sl->coord_system()->is_default());
    print $fh join( "\t", ( $sl->seq_region_name(), $sl->length() ) ), "\n";
  }

  close($fh);
}

1;
