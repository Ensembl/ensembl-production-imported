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

    Bio::EnsEMBL::EGPipeline::PHIBase::PHIAccumulatorWriter

=head1 SYNOPSIS

   
    
=head1 DESCRIPTION

    This is a module for validating all the phi_base candidate entries and grouping them into a super_hash in the accumulator. 

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <https://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <https://www.ensembl.org/Help/Contact>.

=cut


package Bio::EnsEMBL::EGPipeline::PHIBase::PHIAccumulatorWriter;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::Process');
use Data::Dumper;


sub param_defaults {
    return {
        'funnel_branch_code' => 1,
    };
}

sub fetch_input {
    my $self = shift;

    

    my $_species = $self->param_required('_species');    
    unless (defined $_species) {
        die "_species $_species does not exist"; # Will cause job to fail and leave a message in log_message
    }

    my $_dbname = $self->param_required('_dbname');      
    unless (defined $_dbname) {
        die "_dbname $_dbname does not exist for $_species"; # Will cause job to fail and leave a message in log_message
    }

    my $phibase_id = $self->param_required('phi_entry');
    if (! $phibase_id || _remove_spaces($phibase_id) eq '' ) {
         die 'phibase_id is not defined'; # Will cause job to fail and leave a message in log_message
    } 


    my $host_name = $self->param_required('host_name');
    if (! $host_name || _remove_spaces($host_name) eq '' ) {
         die ' host_name is not defined'; # Will cause job to fail and leave a message in log_message
    }

    my $phenotype_name = $self->param_required('phenotype');
    if (! $phenotype_name || _remove_spaces($phenotype_name) eq '' ) {
         die ' phenotype_name is not defined'; # Will cause job to fail and leave a message in log_message
    }

    my $evidence = $self->param_required('_evidence');
    if (! $evidence || _remove_spaces($evidence) eq '' ) {
         die " evidence is not defined for $phibase_id _dbname $_dbname"; # Will cause job to fail and leave a message in log_message
    }

    my $translation_id = $self->param('_translation_id');
    if (! $translation_id || $translation_id eq '' ) {
      die "No translation_id :$translation_id: to write to the accumulator\n$phibase_id :: $_species :: $_dbname\t"; # Will cause job to fail and leave a message in log_message
    } 

}


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;
    
    my $associated_xrefs = $self->_build_associated_xrefs();
    $self->param('associated_xrefs',$associated_xrefs);

}

=head2 write_output

    Description : Implements write_output() interface method of Bio::EnsEMBL::Hive::Process that is used to deal with job's output after the execution.

=cut

sub write_output {  # nothing to write out, but some dataflow to perform:
    my $self = shift @_;
    my $funnel_branch_code = $self->param('funnel_branch_code');
    my $associated_xrefs = $self->param('associated_xrefs');

    my $translation_db_id = $self->param('translation_db_id');
    if (! $translation_db_id || $translation_db_id eq '' ) {
      print "No translation to write to the accumulator\n"; # Will cause job to fail and leave a message in log_message
      return;
    } else {
      my $phi_output_id = $self->_build_output_id_hash($associated_xrefs);
      $self->dataflow_output_id($phi_output_id, $funnel_branch_code);
    }

    
}

=head2 _build_associated_xrefs

    Description : Builds a hash with all the keys that are going to be writen to the DB: 
    'phenotype', 'host', experimental condition', 'DOI', 'Pubmed ref' and 'evidence type'(BLAST/DIRECT MATCH)
    Returns 0 if there is not translation associaed to this entry.

=cut

sub _build_associated_xrefs {  
    my $self = shift @_;

    my $translation_stable_id = $self->param('_translation_id');
    my $phibase_id = $self->param('phi_entry');
    my $phenotype_name = $self->param('phenotype');  
    my $host_name = $self->param('host_name');
    my $host_tax_id = $self->param('host_tax_id');
    my $evidence = $self->param('_evidence');
    my $species = $self->param('_species');
    my $condition_names = $self->param('experiment_condition');
    my $literature_ids = $self->param('litterature_id');
    my $dois = $self->param('DOI');
    my $percent_identity = $self->param('_percent_identity');

    my $associated_xrefs = {};

    if (! $translation_stable_id ) {
         print " ------------------------------\nPHI ACCUMULATOR WRITER :: \ntranslation_stable_id is NOT DEFINED" . 
                                               "\tfor $phibase_id :: $species :: \n ";
         return ;
    } else {
       print "------------------------------\nPHI ACCUMULATOR WRITER \ntranslation: $translation_stable_id " .
                                            "\tfor $phibase_id // :: $species :: \n";
    }
    
    my $translation_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, "Core", "Translation");
    my $db_connection = $translation_adaptor->dbc();
    my $translation = $translation_adaptor->fetch_by_stable_id($translation_stable_id);

    if (! $translation) {
      print ">> Could not find translation $translation_stable_id on Ensemble database: $species \n";
      return;
    }

     #Adding host species to associated xref
    $associated_xrefs->{host} = { 
                    id => $host_tax_id, 
                    label => $host_name };

    #Adding phenotype and condition to associated xref
    $associated_xrefs->{phenotype} = {
                    id => "$phenotype_name",        
                    label => "$phenotype_name" };

    #Adding experimental condition
  
    my $clean_cond_names = lc( _remove_spaces($condition_names) );

    for my $condition_full_name ( split( /;/, $clean_cond_names ) ) {

    if ( ! ($condition_full_name eq "")) {
        print "condition_full_name:". $condition_full_name .":\n";
        push (@{$associated_xrefs->{condition}}, {
                    id => $condition_full_name,             
                    label => $condition_full_name });
      }
    }

    #Adding litterature references
    #my $clean_literature_ids = lc( _remove_spaces($literature_ids) );
    if ( defined $literature_ids && _remove_spaces($literature_ids) ne '') {
      print "Processing literature refs(s) ID:'$literature_ids' \n";
      for my $publication ( split( /;/, $literature_ids ) ) {
        push @{$associated_xrefs->{pubmed}}, $publication;
      }
    } else {
      print "Could not find litterature ID reference for $phibase_id \n";
    }

    #Adding DOI
    if ( defined $dois && $self->_remove_spaces($dois) ne '') {
      print "Processing literature ref(s) DOI:'$dois' \n";
      for my $publication ( split( /;/, $dois ) ) {
        push @{ $associated_xrefs->{doi} }, $publication;
      }
    } else {
      print "Could not find litterature DOI reference for $phibase_id \n";
    }
    
    push @{ $associated_xrefs->{evidence} }, $evidence;
    $associated_xrefs->{translation_stable_id} = $translation_stable_id;
    $associated_xrefs->{percent_identity} = $percent_identity;
    

    $self->param('translation_db_id',$translation->dbID() );
    $db_connection && $db_connection->disconnect_if_idle();

    return $associated_xrefs;
}

=head2 _build_output_id_hash 

    Description: a private method that returns a hash of parameters to be dataflowed to the accumulator:
                    _species
                    translation_db_id
                    phi_entry
                    associated_xrefs

=cut
sub _build_output_id_hash {
    my ($self, $associated_xrefs) = @_;

    my $job_param_hash = {};
  
    # $job_param_hash->{ '_species' } = $self->param('_species');
    # $job_param_hash->{ 'translation_db_id' } = $self->param('translation_db_id');
    # $job_param_hash->{ 'phi_entry' } = $self->param('phi_entry');
    $job_param_hash->{ 'associated_xrefs' } = $associated_xrefs;

    print "{" . $self->param('_species') . "}{" . $self->param('translation_db_id') . "}{" . $self->param('phi_entry') . "}\n";


    return $job_param_hash;
}

=head2 _remove_spaces 

    Description: a private method that returns a string without leading or trailing whitespaces

=cut
sub _remove_spaces {
  my $string = shift;
  $string =~ s/^\s*//;
  $string =~ s/\s*$//;
  return $string;
}

1;
