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

    Bio::EnsEMBL::EGPipeline::PHIBase::ProteinBlaster

=head1 SYNOPSIS

   
    
=head1 DESCRIPTION

    This is a generic RunnableDB module for cleaning all PHI-base associated xrefs from all the dbs in the registry for a particular division. 

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <https://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <https://www.ensembl.org/Help/Contact>.

=cut


package Bio::EnsEMBL::EGPipeline::PHIBase::ProteinBlaster;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::Process');
use Data::Dumper;
use Bio::EnsEMBL::EGPipeline::Xref::BlastSearch;
use Bio::Root::Root;
use Bio::Tools::Run::StandAloneBlastPlus;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::Seq;
use Carp;

sub param_defaults {
    return {

        'flow_branch_code'   => 3,
    };
}

sub fetch_input {
    my $self = shift;
    
    my $phi_entry = $self->param_required('phi_entry');    
    unless (defined $phi_entry) {
        die "phi_entry $phi_entry does not exist"; # Will cause job to fail and leave a message in log_message
    }

    my $annotn_uniprot_acc = $self->param_required('uniprot_acc');
    if (! $annotn_uniprot_acc || _remove_spaces($annotn_uniprot_acc) eq '' ) {
        die 'uniprot accesion is not defined'; # Will cause job to fail and leave a message in log_message
    }

    my $branch_species = $self->param_required('_species');
    if (! $branch_species || _remove_spaces($branch_species) eq '' ) {
        die 'Branch species is not defined'; # Will cause job to fail and leave a message in log_message
    }

    my $db_name = $self->param_required("_dbname");
    if (! $db_name) {
         die "dbname is not defined"; # Will cause job to fail and leave a message in log_message
    }
}


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;
  
    my $annotn_uniprot_acc = $self->param_required('uniprot_acc');
    my $branch_species = $self->param_required('_species');
    my $db_name = $self->param_required("_dbname");
    my $phi_entry = $self->param_required('phi_entry');

    if ($annotn_uniprot_acc =~ /no data/) {
      print "$annotn_uniprot_acc / no uniprot accession for $phi_entry\n";
      return;
    }

    print "-------------------------------- \n$phi_entry\n";


    my $brch_dba = Bio::EnsEMBL::Registry->get_DBAdaptor( $branch_species, "core" );
    if ($brch_dba->dbc()->dbname() =~ /collection/) {
      $brch_dba->{'_is_multispecies'} = 1;
    }
    
    my ($uniprots, $id_type) = $self->get_translation_sequence($annotn_uniprot_acc);

    if ($uniprots->{$annotn_uniprot_acc}->{seq} ne '') {
        print "\t-Got uniprot/uniparc accession.- \n";

        #  Create BLAST DB for this specie if it doesn't already exists or if it needs to be updated
        my $first_char = substr($branch_species, 0, 1); # to avoid too many DBs in a folder we subdivide by qalphabetical order
        my $blast_db_dir = $self->param('blast_db_directory') . "$first_char" ;
        my $blast_database_file = $blast_db_dir . "/" . $branch_species . ".pep.all.fa";
     
        print "BLAST search: creating peptide file if it doesnt exist in $blast_database_file:\n";

        my $pept_file_update_needed = $self->is_update_needed($brch_dba, $branch_species);
     
        if ( ! -f $blast_database_file || $pept_file_update_needed) { # Returns undef if the file doesn't exist, or if the file exists but it's not plain file (e.g. if it's a directory).
          $self->get_peptide_file($brch_dba);    
        } else {
          print "Update not needed, using previously existing DB.\n";
        }

        my $db_connection = $brch_dba->dbc();
        $db_connection && $db_connection->disconnect_if_idle();

        #Create a BlastPlus factory
        my $fac;
    
        eval {
          $fac = Bio::Tools::Run::StandAloneBlastPlus->new( 
            -db_name => $branch_species,
            -db_dir => $blast_db_dir ,
            -db_data => $blast_database_file,
            -create => 1
          );
   
          $fac->make_db();
        };
        $@ =~ /EXCEPTION|WARNING/ and my $e = $@;
        if ( defined $e || ! defined($fac) ) {
          print "**********\nBlastPlus factory could not be created for accession $annotn_uniprot_acc: 
                 Fix this and rerun, delete this entry from the csv file or forgive this job.\n**********\n";
          return 0;
        }

        my $query   = Bio::Seq->new( -display_id => $annotn_uniprot_acc, -seq =>  $uniprots->{$annotn_uniprot_acc}->{seq} );
        my $results = $fac->blastp( -query => $query, -method_args => [ -num_alignments => 1 ]);               
        my $query_length = length($uniprots->{$annotn_uniprot_acc}->{seq});

    
        $self->param('results', $results);
        $self->param('query_length', $query_length);
    }

}

=head2 get_peptide_file
    
    Description: private method that connects to the given species db and fetches all its translatin to a loal file.

=cut

sub get_peptide_file {
    my ($self, $brch_dba) = @_;
    
    my $branch_species = $self->param_required('_species');
    my $first_char = substr($branch_species, 0, 1); # to avoid too many DBs in a folder we subdivide by alphabetical order
    my $blast_db_dir = $self->param('blast_db_directory') . "$first_char" ;
    my $blast_database_file = $blast_db_dir . "/" . $branch_species . ".pep.all.fa";
    my $checksum_file = $blast_db_dir . "/" . $branch_species . ".cks";
    
    open (my $db_file, ">", $blast_database_file) or croak "Could not open ".$blast_database_file;;
    print "update_needed, getting all translations\n";

    my $translation_adaptor=$brch_dba->get_TranslationAdaptor();
    my @translations = @{$translation_adaptor->fetch_all()}; 

    foreach my $translation (@translations) { # store translations into the pep.all.fa file 
      my $stable_id = $translation->stable_id();
      my $sequence = $translation->seq();
      print $db_file ">" . $stable_id . "\n";
      print $db_file $sequence . "\n";
    }
    close $db_file;

    my $checksum = $self->param('checksum');

    open(my $cksfile, '>', $checksum_file); # store fresh checksum value 
    print $cksfile $checksum;
    close $cksfile;
}

=head2 is_update_needed

    Description : private method that decides if a peptide database file for a given species needs to be updated.
    It does so by :
      1) A temporary table with all translations of the species is created in its db.
      2) The checksum of that table is then compared with the checksum value stored with the downloaded peptide file that is being evaluating.
      3) If there is no such file or the checksum don't match the subroutine returns 1, otherwise 'undef';

=cut

sub is_update_needed   {   
    my ($self, $brch_dba, $branch_species) = @_;

    my $pept_file_update_needed = undef;

    if ($self->param('skip_blast_update') == 1 ) { # No update needed
      return 0;
    } 

    my $blast_db_dir = $self->param('blast_db_directory') . substr($branch_species, 0, 1) ;
    my $checksum_file = $blast_db_dir . "/" . $branch_species . ".cks";

    my $dbc=$brch_dba->dbc();
    my $temp_table_name = "ttn_$branch_species";
    if (length($temp_table_name) > 64) {
      $temp_table_name = substr $temp_table_name, 0,64;
    } 
    
    my $sth = $dbc->prepare( "CREATE TEMPORARY TABLE $temp_table_name select  tl.translation_id, tl.modified_date
                                from translation tl join transcript tc using(transcript_id) 
                                join seq_region sr using(seq_region_id) 
                                join coord_system cs using(coord_system_id) 
                                join meta m using(species_id) 
                                where meta_key='species.production_name' and meta_value like '$branch_species';");
    $sth->execute() or
     die "Error creating temporary table for $branch_species: perhaps the DB doesn't have a meta table or the table already exists?\n" .
       "$DBI::err .... $DBI::errstr\n";

    $sth = $dbc->prepare( "checksum TABLE $temp_table_name ;");
    $sth->execute() or
     die "Error obtaining checksum for temporary table $temp_table_name\n" .
       "$DBI::err .... $DBI::errstr\n";

    my $table_name;
    my $checksum;
    $sth->bind_columns(\$table_name,\$checksum);
    $sth->fetch();
    print "Checksum for $branch_species = $checksum\n";
    $self->param('checksum',$checksum);
    
    # Move the delete statement to a bulk loop at the end of db_load_hasher
    $sth = $dbc->prepare( "DROP TABLE $temp_table_name;");
    $sth->execute();

    if ( ! -f $checksum_file ) {
      $pept_file_update_needed = 1;
    } else {
        open my $file, '<', $checksum_file; 
        my $stored_checksum = <$file>; 
        close $file;
        if ($stored_checksum ne $checksum ) {
            print "$checksum is different than expected checksum ($stored_checksum). Blast peptide fasta file update needed\n";
            $pept_file_update_needed = 1;
        }
    }
    
    return $pept_file_update_needed;
}

=head2 write_output

    Description : Implements write_output() interface method of Bio::EnsEMBL::Hive::Process that is used to deal with job's output after the execution.
    param('flow_branch_code'): defines the branch where the fan of jobs is created (1 by default).
    If there is no translation_stable_id supplied, no job will be created

=cut

sub write_output {  
    my $self = shift @_;
    
    my $branch_species = $self->param('_species');
    my $results = $self->param('results');
    
    if (!$results) {
      print "**********\n$branch_species doesn't have results\n";
      return;
    }
    


    my $hit = $results->next_hit();

    
    my $query_length = $self->param_required('query_length');
    my $current_identity_thr = $self->param_required('min_blast_identity');
    my $current_query_length_ratio = $self->param_required('current_query_length_ratio');
    my $flow_branch_code = $self->param('flow_branch_code');

    #check for hit results before going any further
    if ( $results && $hit ) {
      my $translation_stable_id = $hit->name(); 
      if (!$translation_stable_id || $translation_stable_id eq '') {
        return;
      }

      my $hit_length = $hit->length();      
      my $identity = $hit->hsp->percent_identity();

      my $max_length = ($hit_length, $query_length)[$hit_length < $query_length];
      my $min_length = ($hit_length, $query_length)[$hit_length > $query_length];
      my $query_length_covered = $min_length / $max_length ; # inforce 0 < ratio < 1
      
      if ( ($identity >= $current_identity_thr ) && ($query_length_covered >= $current_query_length_ratio)  ) {
        my $translation_params = $self->_build_output_hash($translation_stable_id, $identity, $hit_length );
        $self->dataflow_output_id($translation_params, $flow_branch_code);
      } else {
        return
      }

      print "IDENTITY : $identity ; query-target RATIO: $query_length_covered  \n";
    }

}

=head2 get_translation_sequence 

    Description: a private method that returns the protein sequence corresponding to a given accession.
    Insure that uniprot/uniparc queries only one protein sequence and no more. 

=cut

sub get_translation_sequence {
  my ($self, $annotn_uniprot_acc) = @_;

  my $id_type;
  my $uniprots = {};
  my $up;

  eval {
    $up = $self->get_uniprot_seq($annotn_uniprot_acc);
    $uniprots->{$annotn_uniprot_acc} = $up;
    $id_type = 'uniprot';
  };

  if($@ || $uniprots->{$annotn_uniprot_acc}->{seq} eq '') {
    eval {
      $up = $self->get_uniparc_seq($annotn_uniprot_acc);
      $uniprots->{$annotn_uniprot_acc} = $up;
      $id_type = 'uniparc';
    } 
  }

  
  $self->param('id_type', $id_type);
  # Insure that uniprot/uniparc queries only one protein sequence and no more. 
  my $count_seqs = $up->{seq}   =~ tr/>//; 
  if ($count_seqs > 0) {
    my $annotn_phi_entry = $self->param_required('phi_entry');
    die  " *-*-*- ERROR : skipping $annotn_phi_entry  *-*-* Uniprot/Uniparc query returns more than one sequence\n";
  }

  return ($uniprots);
}


=head2 get_uniprot_seq 

    Description: a private method that returns a protein sequence (or more!) from a accession query to UNIPROT

=cut

sub get_uniprot_seq {
  my ($self, $acc) = @_;
  
  my $search = Bio::EnsEMBL::EGPipeline::Xref::BlastSearch->new();
  my $seq = $search->get('https://www.uniprot.org/uniprot/'.$acc.'.fasta');
  $seq =~ s/^>\S+\s+([^\n]+)\n//;
  my $des = $1;
  $seq =~ tr/\n//d;
  return {seq=>$seq,des=>$des};
}

=head2 get_uniparc_seq 

    Description: a private method that returns a protein sequence (or more!) from a accession query to UNIPARC

=cut
sub get_uniparc_seq {
  my ($self, $acc) = @_;
  
  my $search =  Bio::EnsEMBL::EGPipeline::Xref::BlastSearch->new();
  my $seq = $search->get('https://www.uniprot.org/uniparc/?query=' . $acc . '&columns=sequence&format=fasta');
  $seq =~ s/^>\S+\s+([^\n]+)\n//;
  my $des = $1;
  $seq =~ tr/\n//d;
  return {seq=>$seq,des=>$des};
}

=head2 _build_output_hash 

    Description: a private method that returns a hash of parameters to add to the input plus() 

=cut

sub _build_output_hash {
    my ($self, $translation_stable_id, $percent_identity) = @_;
    # line_fields : All found in SubtaxonsFinder->_substitute_rows  + 
    #                       -------- UPDATED FIELDS -----------    #                      

    #                      '_id_type' # 18 Protein DB (UNIPROT or UNIPARC) used to get the sequence blast_match
    #                      '_translation_id' # 19 translation stable _id
    #                      '_percent_identity' # 20 percentage identity of the match
    
    my @hashes = ();
    my $job_param_hash = {};

    my $id_type = $self->param_required("id_type");
    if (! $id_type) {
         die "id_type is not defined"; # Will cause job to fail and leave a message in log_message
    }

    print "Blasted translation:" . $translation_stable_id . ": found in $id_type, with identity: $percent_identity \n";

    $job_param_hash->{ '_id_type' } = $id_type;
    $job_param_hash->{ '_translation_id' } = $translation_stable_id;
    $job_param_hash->{ '_percent_identity' } = $percent_identity;

    push @hashes, $job_param_hash;
    return \@hashes;
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
