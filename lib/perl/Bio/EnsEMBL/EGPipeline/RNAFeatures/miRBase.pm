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

package Bio::EnsEMBL::EGPipeline::RNAFeatures::miRBase;

use strict;
use warnings;

use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Path::Tiny qw(path);

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub run {
  my ($self) = @_;
  my $species    = $self->param_required('species');
  my $logic_name = $self->param_required('logic_name');
  my $db_name    = $self->param_required('db_name');
  my $files      = $self->param_required('files');
  
  if (exists $$files{$species}) {
    my $dba  = $self->core_dba();
    my $sa   = $dba->get_adaptor('Slice');
    my $dafa = $dba->get_adaptor('DnaAlignFeature');
    my $aa   = $dba->get_adaptor('Analysis');
    
    my $analysis       = $aa->fetch_by_logic_name($logic_name);
    my $external_db_id = $self->fetch_external_db_id($db_name);
    
    $self->parse_mirbase($$files{$species}, $sa, $dafa, $analysis, $external_db_id);
  }
}

sub parse_mirbase {
  my ($self, $file, $sa, $dafa, $analysis, $external_db_id) = @_;
  
  my $path = path($file);
  my $mirbase = $path->slurp;
  
  my @lines = $mirbase =~ /(.*\tmiRNA_primary_transcript\t.*)/gm;
  
  foreach my $line (@lines) {
    chomp $line;
    my (
      $seqname,
      undef,
      undef,
      $start,
      $end,
      undef,
      $orientation,
      undef,
      $attributes
    ) = split(/\t/, $line);
  
    $seqname =~ s/^chr//i;
    my $slice = $sa->fetch_by_region('toplevel', $seqname);
    
    if (defined $slice) {
      my $strand = $orientation eq '-' ? -1 : 1;
      
      my ($accession, $rna_name) = $attributes =~ /ID=([^;]+).*Name=([^;]+)/;
      
      my $length = $end - $start + 1; 
      my $cigar = $length.'M';

      my $biotype = "pre_miRNA";
      
      my @attribs;
      push @attribs, $self->create_attrib('rna_gene_biotype', $biotype);
      push @attribs, $self->create_attrib('rfam_accession',   $accession);
      # I know that says 'rfam_accession', but it works for any accession;
      # the attrib_type will be changed in due course.
      
      my $feature = Bio::EnsEMBL::DnaDnaAlignFeature->new(
        -slice          => $slice,
        -start          => $start,
        -end            => $end,
        -strand         => $strand,
        -hstart         => 1,
        -hend           => $length,
        -hstrand        => 1,
        -hseqname       => $rna_name,
        -cigar_string   => $cigar,
        -external_db_id => $external_db_id,
      );
      $feature->add_Attributes(@attribs);
      $feature->analysis($analysis);
      
      $dafa->store($feature);
    } else {
      $self->throw("miRBase data in $file annotated on unrecognised sequence name");
    }
  }
}

sub create_attrib {
  my ($self, $key, $value) = @_;
  
  my $attrib = Bio::EnsEMBL::Attribute->new(
    -CODE  => $key,
    -VALUE => $value
  );
  
  return $attrib;
}

1;
