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

    Bio::EnsEMBL::EGPipeline::PHIBase::TranslationFinder

=head1 SYNOPSIS

    

=head1 DESCRIPTION

    This runnable processes the PHI_base entry for a particular genome (species taxonomic  node), 
    and tries to find the translation reported by the input file.

    It first tries to do so by finding a direct match in Ensembl of the accession number, gene name or gene locci.
    If a direct match is not found the accession number is used toretrieve the sequence from UNIPROT or UNIPARC.
  
    The sequence is then blasted against the protein db of the particular species.


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <https://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <https://www.ensembl.org/Help/Contact>.

=cut


package Bio::EnsEMBL::EGPipeline::PHIBase::TranslationFinder;

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::EGPipeline::Xref::BlastSearch;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::LookUp::LocalLookUp;
use Bio::EnsEMBL::DBSQL::TaxonomyNodeAdaptor;


use base ('Bio::EnsEMBL::Hive::Process');
use Data::Dumper;
use Time::HiRes qw( time );

use Scalar::Util qw(reftype);

sub param_defaults {
        return {
        'accumulator_branch_code' => 1,
        'blast_branch_code' => 3,
        'skip_blast_branch_code' => 4,
        'phibase_xrefs' => undef,
    };
}


sub fetch_input {
    my $self = shift;

    my $annotn_tax_id = $self->param_required("species_tax_id");
    if (! $annotn_tax_id) {
         die 'annotn_tax_id is not defined'; # Will cause job to fail and leave a message in log_message
    } 

    my $specific_species_name = $self->param_required("_species");
    if (! $specific_species_name) {
         die "specific_species_name ('_species')is not defined"; # Will cause job to fail and leave a message in log_message
    }
    
    my $dbname = $self->param_required("_dbname");
    if (! $dbname) {
         die "dbname is not defined"; # Will cause job to fail and leave a message in log_message
    }

    my $skip_blast_update = $self->param('skip_blast_update');

}


sub run {
    my $self = shift;

    my $specific_species_name = $self->param_required("_species");
    my $translation;
    my $id_type;
    my $brch_dba;
    ($translation, $id_type, $brch_dba) = $self->_find_translation();
    
    $self->param('translation', $translation);
    $self->param('id_type', $id_type);
    $self->param('brch_dba', $brch_dba);
  
}

sub write_output {
    my $self = shift;

    my $blast_branch_code = $self->param('blast_branch_code');
    my $accumulator_branch_code = $self->param('accumulator_branch_code');
    my $skip_blast_branch_code = $self->param('skip_blast_branch_code');
    my $skip_blast_update = $self->param('skip_blast_update');
    my $translation = $self->param('translation');

    my $id_type = $self->param('id_type');
    my $brch_dba = $self->param('brch_dba');

    my $translation_found_params ;
    my $translation_unknown_params_BLAST ;
    my $translation_unknown_params_SKIP_UPDATE ;

    if($translation) {
      print "Translation_found:" . $translation->stable_id . ":\n";
      $translation_found_params = $self->_build_output_hash($translation, $id_type, $brch_dba);
      $self->dataflow_output_id($translation_found_params, $accumulator_branch_code);
      
    } elsif (!$skip_blast_update) {
      print "Blast translation_unknw:$translation:\n";
      $translation_unknown_params_BLAST = $self->_build_output_hash($translation, $id_type, $brch_dba);
      $self->dataflow_output_id($translation_unknown_params_BLAST, $blast_branch_code);
    }  else {
      print "Blast translation_unknw:$translation: (but skipping blastDB update)\n";
      $translation_unknown_params_SKIP_UPDATE = $self->_build_output_hash($translation, $id_type, $brch_dba);
      $self->dataflow_output_id($translation_unknown_params_SKIP_UPDATE, $skip_blast_branch_code);
    }
    
}



=head2 _get_lookup
    
    Description: a private method that loads the registry and the lookup taxonomy DB.

=cut

sub _get_lookup {
    my $self = shift @_;
    
    my $lookup = Bio::EnsEMBL::LookUp::LocalLookUp->new( -SKIP_CONTIGS => 1,
                                                         -NO_CACHE     => 1 ); # Used to check all leafs sub_specie/strains under a taxonomy_id (specie)

    my $db_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor( "multi", "taxonomy" );
    my $tax_adaptor = Bio::EnsEMBL::DBSQL::TaxonomyNodeAdaptor->new($db_adaptor  );
    $lookup->taxonomy_adaptor($tax_adaptor);

    return ($lookup);
}

=head2 _find_translation
    
    Description: a private method that tries to find a translation in the Ensemble DB by either protein_accession, locus or by gene_type.

=cut

sub _find_translation { 
  my $self = shift;
  
  my $specific_species_name = $self->param_required("_species");
  my $specific_db_name = $self->param_required("_dbname");
  my $phi_entry = $self->param_required('phi_entry');   

  print "-----------------------------\n$phi_entry ::  $specific_species_name in $specific_db_name  ---------------------------------\n";

  my $lookup = $self->_get_lookup();
  
  #build dbadaptor from specific_db_name and specific_species_name. This is messy but I didn't find another way of passing the dba along the analysis.
  #If you find a better way of doing please improve this
  my $dba = Bio::EnsEMBL::Registry->get_DBAdaptor( $specific_species_name, "core" );

  my $uniprot_acc = $self->param_required('uniprot_acc');
  my $locus_tag = $self->param_required('locus');
  my $gene_name = $self->param_required('gene_name');
  print "-----------------------------  accession $uniprot_acc /locus $locus_tag /gene_name $gene_name ---------------------------------\n";

  my $translation = '';
  my $identifier_type = '';

  my $dbentry_adaptor = $dba->get_adaptor("DBEntry");
  my @transcripts_ids = $dbentry_adaptor->list_transcript_ids_by_extids($uniprot_acc); 

  print " nb transcript ids by uniprot:" . scalar(@transcripts_ids) . "\n";
  my @gene_ids = $dbentry_adaptor->list_gene_ids_by_extids($uniprot_acc); # List of gene_id by an external identifier accession that is
                                                                          # linked to  any of the genes transcripts, translations or the gene itself
  print " nb gene_ids by uniprot:" . scalar(@gene_ids) . "\n";                                                                          

  if ( scalar(@gene_ids) == 0 ) {
    @gene_ids = $dbentry_adaptor->list_gene_ids_by_extids($locus_tag);   # List of gene_id by an external identifier locus that is
                                                                          # linked to  any of the genes transcripts, translations or the gene itself
    $identifier_type = 'locus';                                               

    print " nb gene_ids by locus:" . scalar(@gene_ids) . "\n";
  }
  if ( scalar(@gene_ids) == 0 ) {
    @gene_ids = $dbentry_adaptor->list_gene_ids_by_extids($gene_name);   # List of gene_id by an external identifier gene name that is
                                                                          # linked to  any of the genes transcripts, translations or the gene itself
    $identifier_type = 'name';
    print " nb gene_ids by name:" . scalar(@gene_ids) . "\n";
  }

  my $translation_adaptor = $dba->get_adaptor("Translation");
  my $transcript_adaptor  = $dba->get_adaptor("Transcript");
  my $gene_adaptor        = $dba->get_adaptor("Gene");
  my $transcript;
  my $gene;
  
  if ( scalar(@transcripts_ids) >= 1 ) {
    
    my $transcript_id = $transcripts_ids[0];
    $transcript = $transcript_adaptor->fetch_by_dbID($transcript_id);
    $translation = $translation_adaptor->fetch_by_Transcript($transcript);
    $identifier_type = 'accession';
  } elsif ( scalar(@gene_ids) >= 1 ) {
      $gene = $gene_adaptor->fetch_by_dbID( $gene_ids[0] );
      my @transcripts = @{ $transcript_adaptor->fetch_all_by_Gene($gene) };
      $transcript = $transcripts[0];
      $translation = $translation_adaptor->fetch_by_Transcript($transcript);
  } else {
      print "NO TRANSCRIPTS NO GENES FOR THIS:: uniprot_acc: $uniprot_acc:, locus_tag: $locus_tag:, gene_name: $gene_name:, identifier_type: $identifier_type;\n";
  }

  return $translation, $identifier_type, $dba;
} 


=head2 _accumulate_phibase_xrefs 

    Description: a private method that, in those cases where the translation is found via a DIRECT_MATCH'
    builds a hash containing the phibase_xrefs and pass it to the accumulator.

=cut

sub _accumulate_phibase_xrefs {
    my $self = shift;

    my $phibase_id = $self->param_required('phi_entry');
    if (! $phibase_id || remove_spaces($phibase_id) eq '' ) {
         die 'phibase_id is not defined'; # Will cause job to fail and leave a message in log_message
    } 

    my $acc = $self->param_required('uniprot_acc');
    if (! $acc || remove_spaces($acc) eq '' ) {
         die 'uniprot accesion is not defined'; # Will cause job to fail and leave a message in log_message
    } 

    my $host_tax_id = $self->param_required('host_tax_id');
    if (! $host_tax_id || remove_spaces($host_tax_id) eq '' ) {
         die ' host_tax_id is not defined'; # Will cause job to fail and leave a message in log_message
    }

    my $host_name = $self->param_required('host_name');
    if (! $host_name || remove_spaces($host_name) eq '' ) {
         die ' host_name is not defined'; # Will cause job to fail and leave a message in log_message
    }

    my $phenotype_name = $self->param_required('phenotype');
    if (! $phenotype_name || remove_spaces($phenotype_name) eq '' ) {
         die ' phenotype_name is not defined'; # Will cause job to fail and leave a message in log_message
    }

    my $evidence = $self->param_required('_evidence');
    if (! $evidence || remove_spaces($evidence) eq '' ) {
         die ' evidence is not defined'; # Will cause job to fail and leave a message in log_message
    }

    # nested hash of phibase xrefs by species-translation_id-phiId
    my $phibase_xrefs = {};

}

=head2 _build_output_hash 

    Description: a private method that returns a hash of parameters to add to the input plus() 

=cut

sub _build_output_hash {
    my ($self, $translation, $id_type, $dba) = @_;
    my @hashes = ();

    # line_fields : All found in SubtaxonsFinder->_substitute_rows  + 
    #                       -------- NEW FIELDS -----------    #                      
    #                      '_evidence' # 17 is this a direct_match or a blast_match
    #                      '_id_type' # 18 How the match was found (if DIRECT_MATCH locus_tag, gene name or accession), or Protein DB used to get the sequence blast_match
    #                      '_translation_id' # 19 translation stable _id
    #                      '_percent_identity' # 20 percentage identity of the match
    #

    my $job_param_hash = {};
    my $evidence = 'SEQUENCE_MATCH';

    if ($translation) {
      print "\t\t--- Found gene on ensembl: ";
      my $gene_adaptor = $dba->get_adaptor("Gene");
      my $db_connection = $gene_adaptor->dbc();
      my $gene = $gene_adaptor->fetch_by_translation_stable_id($translation->stable_id);
      
      $db_connection && $db_connection->disconnect_if_idle();
      print "using $id_type\n\t\t" . $gene->stable_id . "\t (should have a translation):" . $translation->stable_id . ":\n";
      
      $evidence = 'DIRECT';
      $job_param_hash->{ '_id_type' } = $id_type;
      $job_param_hash->{ '_translation_id' } = $translation->stable_id;
      $job_param_hash->{ '_percent_identity' } = 100;

    } else {
        print "\t\tDidn't find direct match on Ensembl \t id_type :$id_type\n";

        $job_param_hash->{ '_id_type' } = undef;
        $job_param_hash->{ '_translation_id' } = undef;
        $job_param_hash->{ '_percent_identity' } = undef;
      }

    $self->param('evidence', $evidence);
    $job_param_hash->{ '_evidence' } = $evidence ;

    push @hashes, $job_param_hash;
    
    return \@hashes;
  
}

1;
