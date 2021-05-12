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

    Bio::EnsEMBL::EGPipeline::PHIBase::SubtaxonsFinder

=head1 SYNOPSIS


=head1 DESCRIPTION

    This runnable processes a particular PHI-base entry from the inputfile, and from its host species
    queries the taxonomy db to find its taxonomic children. It fans one job per subspecies returned.
    In the event it finds too many branches, it returns a restricted selection of sub-branches. 
    A parameter (MAX_SUB_TAX_DBAS) defines the maximum number of taxonomic children to accept. 

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <https://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <https://www.ensembl.org/Help/Contact>.

=cut


package Bio::EnsEMBL::EGPipeline::PHIBase::SubtaxonsFinder;

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::LookUp::LocalLookUp;
use Bio::EnsEMBL::DBSQL::TaxonomyNodeAdaptor;

use base ('Bio::EnsEMBL::Hive::Process');
use Data::Dumper;
use Time::HiRes qw( time );

use Scalar::Util qw(reftype);

sub param_defaults {
        return {
        'MAX_SUB_TAX_DBAS' => 30,
        'fan_branch_code'  => 2,
        'has_dba'          => undef,
    };
}

sub pre_cleanup {
    # This is where the runnable would take care of getting things
    # (typically files or database connections) in order before starting to run.
    # This is often useful if the jobs sometimes need to be retried before
    # completing successfully
}

sub fetch_input {
    my $self = shift;
    # read_input is typically used to validate input, test db connections are
    # working, etc. 
    my $annotn_tax_id = $self->param_required("species_tax_id");
    if (! $annotn_tax_id) {
         die 'annotn_tax_id is not defined'; # Will cause job to fail and leave a message in log_message
    }
}


sub run {
    my $self = shift;
    
    my $annotn_tax_id = $self->param_required("species_tax_id");
    my $phi_entry = $self->param_required('phi_entry');    
    my $pathogen_name = $self->param_required('pathogen_name');  
    print "--------------------------- phi_entry: $phi_entry  :: $pathogen_name with tax_id $annotn_tax_id\n";

    
    my $output_subtaxons =  $self->_get_subtaxons_dbas($annotn_tax_id);
    
    $self->param('subtaxons_dbas', $output_subtaxons);
}

=head2 _get_subtaxons_dbas
    
    Description: a private method that returns an arrayref of all dbadaptors of a supplied taxonomy_id ( = with its taxonomic children if any).
    In the eventuality of too many taxonomic children, the methood returns a maximum number of subtax dbas, limit defined by MAX_SUB_TAX_DBAS

=cut
sub _get_subtaxons_dbas {
  my ($self, $annotn_tax_id) = @_;
  my $lookup = $self->_get_lookup();
  my $db_connection = $lookup->taxonomy_adaptor->dbc();

  my $branch_dbas = $lookup->get_all_by_taxon_branch($annotn_tax_id);
  my $nb_taxon_branch_dbas = scalar(@{$branch_dbas});
  my $MAX_SUB_TAX_DBAS = $self->param_required("MAX_SUB_TAX_DBAS");

  print "nb tax dbs: $nb_taxon_branch_dbas \n";

    if ($nb_taxon_branch_dbas == 0) {
      print "\tNo dbs for this specie tax_id: $annotn_tax_id\n";
  } elsif ( $nb_taxon_branch_dbas > $MAX_SUB_TAX_DBAS ) {# if too many branch_dbas, limit the selection to MAX_SUB_TAX_DBAS
      print "\tToo many dbs for this specie (tax_id:$annotn_tax_id):  $nb_taxon_branch_dbas\n";
      $branch_dbas = $self->_limit_branch_dbas(@{$branch_dbas}); 
  } 

  #close de dba connection if it is still open (avoids too many connections error)
  $db_connection && $db_connection->disconnect_if_idle();

  return $branch_dbas;
}

=head2 _limit_branch_dbas
    
    Description: a private method that limits a supplied arrayref of dbas to a predefined number (MAX_SUB_TAX_DBAS).
    The priority to select among a list of dbas is:
    1- If a subtaxon_id is reported in the input csv file, include it, it's the one experimentally verified.
    2- If a core_db is among the list of dbas, include it. So far only core DBs are searchable with Biomart.
    3- Keep adding among the remaining dbas until reaching MAX_SUB_TAX_DBAS limit.

=cut

sub _limit_branch_dbas {
  my ($self, @branch_dbas) = @_;
  my $new_branch_dbas;
  my $nb_selected_dbas = 0;
  my $nb_taxon_branch_dbas = scalar(@branch_dbas);
  my $MAX_SUB_TAX_DBAS = $self->param_required("MAX_SUB_TAX_DBAS");
  my $lookup = $self->param('lookup');
  my $species_tax_id = $self->param('species_tax_id'); #species taxonomic node level
  my $subtaxon_id = $self->param('subtaxon_id'); # strain taxa node level
  my $reported_dba;
  my %selected_species;

  # if strain level tax_id is provided by input file, make sure its dba(s) is returned.
  if ($subtaxon_id && $subtaxon_id ne '' && $subtaxon_id ne 'no data found' && $subtaxon_id != $species_tax_id) {
    $reported_dba = $lookup->get_all_by_taxon_id($subtaxon_id);

    foreach my $rdba (@$reported_dba) { # might return more than one dba      
      my $species_name = $rdba->{'_species'};
      if ($species_name ){ # Do not accept an empty node
        push @{$new_branch_dbas}, $rdba;
        $selected_species{$species_name} = 1;
        $nb_selected_dbas ++;
      }
    }
  }

  # Go through all dbs and select all coreDBs (in exceptional cases there might be more than one coreDB but most often there will be only one)
  for (my $i = 0; $i < $nb_taxon_branch_dbas; $i++) {
    my $brch_dba = $branch_dbas[$i]; 
    my $species_name = $brch_dba->{'_species'};
    if ( ($brch_dba->dbc()->dbname() !~ /collection/)  &&  (!$selected_species{$species_name}) ) { 
      push @{$new_branch_dbas}, $brch_dba;
      $nb_selected_dbas ++;
    } 
  }
  
  # Now add a limited number of collection DBs (minus the already selected)
  for (my $i=0; $i < ($MAX_SUB_TAX_DBAS-$nb_selected_dbas) && $i < $nb_taxon_branch_dbas ; $i++) {
    my $brch_dba = $branch_dbas[$i];
    my $species_name = $brch_dba->{'_species'}; 
    if (($brch_dba->dbc()->dbname() =~ /collection/)   &&  (!$selected_species{$species_name})) {
      push @{$new_branch_dbas}, $brch_dba;
    }
  }  
 
  return $new_branch_dbas;
}

sub write_output {
    my ($self, @branch_dbas) = @_;

    my $fan_branch_code = $self->param('fan_branch_code');
    my $subtax_entries = $self->_build_output_hash();
    
    # "fan out" into fan_branch_code:
    $self->dataflow_output_id($subtax_entries, $fan_branch_code);
    
}

=head2 _build_output_hash 

    Description: a private method that returns a hash of parameters for all subjobs to fan. 
    Each job will have its own corresponding dba + the current job ($self) input fanned parameters from phiFileReader

=cut

sub _build_output_hash {
    my $self = shift;
    my @hashes = ();
    my $subtax_dbas = $self->param_required("subtaxons_dbas");
    # line_fields : All found in phiFileReader->_substitute_rows  + 

    #                     ---- NEW FIELDS ----
    #                    # '_species' # 15, specific species name of the subtaxon branch
    #                    # '_dbname' # 17 subtaxon db name
    foreach my $dba (@$subtax_dbas) {
        my $job_param_hash = {};
        $job_param_hash->{ '_species' } = $dba->{ '_species' };
        $job_param_hash->{ '_dbname' } = $dba->dbc()->{'_dbname'};
        push @hashes, $job_param_hash;
    }

    return \@hashes;
}

sub post_cleanup {

}

=head2 _get_lookup
    
    Description: a private method that loads the registry and the lookup taxonomy DB.

=cut

sub _get_lookup {
    my $self = shift @_;
    
    my $lookup;
    my $wait_and_retry = 1;
    my $max_retry = 20;
    
    # This wrap around the lookup avoids to a certain extent "too many connections" errors, by sleeping and retrying 15 seconds later.
    while ($wait_and_retry && $max_retry) {
      eval {
        $lookup = Bio::EnsEMBL::LookUp::LocalLookUp->new( -SKIP_CONTIGS => 1,
                                                          -NO_CACHE     => 1 ); # Used to check all leafs sub_specie/strains under a taxonomy_id (specie)
        $wait_and_retry = 0;
      };
      if ($@) {
        sleep 15;
        $max_retry = $max_retry - 1;
      }
    }
    if (!$lookup) { die "The LocalLookUp failed to be created\n"};
    my $db_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor( "multi", "taxonomy" );
    my $tax_adaptor = Bio::EnsEMBL::DBSQL::TaxonomyNodeAdaptor->new($db_adaptor  );
    $lookup->taxonomy_adaptor($tax_adaptor);
    $self->param('lookup',$lookup);

    return ($lookup);
}


1;
