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

package Bio::EnsEMBL::EGPipeline::FileDump::GTFDumper;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::FileDump::BaseDumper');

use Bio::EnsEMBL::Utils::IO::GTFSerializer;

use Path::Tiny qw(path);

sub param_defaults {
  my ($self) = @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'data_type'        => 'basefeatures',
    'file_type'        => 'gtf',
    'logic_name'       => [],
  };
}

sub run {
  my ($self) = @_;
  my $species     = $self->param_required('species');
  my $db_type     = $self->param_required('db_type');
  my $out_file    = $self->param_required('out_file');
  my $logic_names = $self->param_required('logic_name');
  
  my $reg = 'Bio::EnsEMBL::Registry';
  
  my $sa = $reg->get_adaptor($species, $db_type, 'Slice');
  my $slices = $sa->fetch_all('toplevel');
  
  open(my $out_fh, '>', $out_file) or $self->throw("Cannot open file $out_file: $!");
  my $serializer = Bio::EnsEMBL::Utils::IO::GTFSerializer->new($out_fh);
  $serializer->print_main_header($self->get_DBAdaptor($db_type));
  
  foreach my $slice (@$slices) {
    my @genes;
    
    if (scalar(@$logic_names) == 0) {
      @genes = @{$slice->get_all_Genes(undef, undef, 1)};
    } else {
      foreach my $logic_name (@$logic_names) {
        push @genes, @{$slice->get_all_Genes($logic_name, undef, 1)};
      }
    }
    
    foreach my $gene (@genes) {
      $serializer->print_Gene($gene);
    }
  }
  
  close($out_fh);
  
  $self->remove_attributes($out_file);
}

sub remove_attributes {
  my ($self, $out_file) = @_;
  
  my @unwanted_atts = qw(
    \w+_biotype
    \w+_name
    \w+_source
    \w+_version
    exon_id
    protein_id
  );
  my $att_regex = "(" . join("|", @unwanted_atts) . ")";
  
  my $tmp_file = $out_file . ".tmp";
  open my $in, "<", $out_file or die "Can't read $out_file\n";
  open my $tmp, ">", $tmp_file or die "Can't read $tmp_file\n";
  for my $line (readline $in) {
    $line =~ s/\s*$att_regex\s*"[^"]+";//g;
    print $tmp $line;
  }
  close $in;
  close $tmp;
  
  rename $tmp_file, $out_file;
}

1;
