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

=head1 NAME

Bio::EnsEMBL::EGPipeline::Xref::LoadGPR

=head1 DESCRIPTION

Add UniParc xrefs to a core database, based on checksums.

=head1 Author

Dan Bolser, Dan Staines and James Allen

=cut

package Bio::EnsEMBL::EGPipeline::Xref::LoadGPR;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::Xref::LoadXref');

use Bio::EnsEMBL::DBSQL::DBAdaptor;

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    'logic_name'  => 'gramene_plant_reactome',
    'uppercase_gene_id' => 0,
  };
}

sub run {
  my ($self) = @_;

  my $db_type      = $self->param_required('db_type');
  my $logic_name   = $self->param_required('logic_name');
  my $external_dbs = $self->param_required('external_dbs');

  my $dba = $self->get_DBAdaptor($db_type);
  my $aa  = $dba->get_adaptor('Analysis');

  my $analysis = $aa->fetch_by_logic_name($logic_name);

  ## What does this do? Reset the version number?
  for my $external_db (keys %$external_dbs){
      $self->external_db_reset($dba, $external_db);
  }

  $self->add_xrefs($dba, $analysis, $external_dbs);

  ## What does this do? Update the version number?
  for my $external_db (keys %$external_dbs){
      $self->external_db_update($dba, $external_db);
  }
}

sub add_xrefs {
  my ($self, $dba, $analysis, $external_dbs) = @_;

  my $reac_file = $self->param_required('xref_reac_file');
  my $path_file = $self->param_required('xref_path_file');

  open RXREF, '<', $reac_file
      or die "Oh for the love of crack!\n";

  open PXREF, '<', $path_file
      or die "Oh for the love of crack!\n";

  my %gene_id_to_reac;
  my %gene_id_to_path;

  while (<RXREF>){
      chomp;

      my ($gene_id,                 # Col 1: Gene_ID
          $reac_id,                 # Col 2: PR_Stable_ID
          undef,                    # Col 3: Direct URL to reaction in browser
          $reac_name,               # Col 4: Reaction Name
          $evidence_code,           # Col 5: Evidence Code
          $species_scientific_name, # Col 6: Species (name)
         ) = split/\t/;

      push @{ $gene_id_to_reac{$gene_id} }, [$reac_id, $reac_name];
  }

  while (<PXREF>){
      chomp;

      my ($gene_id,                 # Col 1: Gene_ID
          $path_id,                 # Col 2: Pathway_ID
          undef,                    # Col 3: Direct URL to reaction in browser
          $path_name,               # Col 4: Pathway Name
          $evidence_code,           # Col 5: Evidence Code
          $species_scientific_name, # Col 6: Species (name)
         ) = split/\t/;

      push @{ $gene_id_to_path{$gene_id} }, [$path_id, $path_name];
  }

  ## Old style (single xref file with both reactions and pathways)

  # while (<XREF>){
  #     chomp;

  #     my ($gene_id,                 # Col 1: Gene_ID
  #         $reac_id,                 # Col 2: PR_Stable_ID
  #         undef,                    # Col 3: Gene_ID [with cellular localization]
  #         $path_id,                 # Col 4: Pathway_ID
  #         undef,                    # Col 5: Direct URL to pathway in browser
  #         $path_name,               # Col 6: Pathway Name
  #         $evidence_code,           # Col 7: Evidence Code
  #         $species_scientific_name, # Col 8: Species (name)
  #        ) = split/\t/;

  #     push @{ $gene_id_to_reac{$gene_id} }, $reac_id;
  #     push @{ $gene_id_to_path{$gene_id} }, [$path_id, $path_name];
  # }

  my $ga   = $dba->get_adaptor('Gene');
  my $dbea = $dba->get_adaptor('DBEntry');

  my $genes = $ga->fetch_all();

  foreach my $gene (@$genes) {
      my $gene_stable_id = $gene->stable_id;
      $gene_stable_id = uc($gene_stable_id) if ($self->param('uppercase_gene_id'));

      my $reac_aref = $gene_id_to_reac{$gene_stable_id};
      my $path_aref = $gene_id_to_path{$gene_stable_id};

      if (defined $reac_aref) {

          for my $reac (@$reac_aref){
              my $xref = $self->
                  add_reac_xref($reac, $analysis, $external_dbs->{'reac'});

              $dbea->store($xref, $gene->dbID(), 'Gene');
          }

          for my $path (@$path_aref){
              my $xref = $self->
                  add_path_xref($path, $analysis, $external_dbs->{'path'});

              $dbea->store($xref, $gene->dbID(), 'Gene');
          }
      }
  }
}

sub add_reac_xref {
  my ($self, $reac, $analysis, $external_db) = @_;

  my $xref = Bio::EnsEMBL::DBEntry->
      new( -PRIMARY_ID => $reac->[0],
           -DISPLAY_ID => $reac->[1],
           -DBNAME     => $external_db,
           -INFO_TYPE  => 'DIRECT',
      );

  $xref->analysis($analysis);

  return $xref;
}

sub add_path_xref {
  my ($self, $path, $analysis, $external_db) = @_;

  my $xref = Bio::EnsEMBL::DBEntry->
      new( -PRIMARY_ID => $path->[0],
           -DISPLAY_ID => $path->[1],
           -DBNAME     => $external_db,
           -INFO_TYPE  => 'DIRECT',
      );

  $xref->analysis($analysis);

  return $xref;
}

1;
