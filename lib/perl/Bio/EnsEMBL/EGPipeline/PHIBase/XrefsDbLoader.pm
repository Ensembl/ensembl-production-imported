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

    Bio::EnsEMBL::EGPipeline::PHIBase::XrefsDbLoader

=head1 SYNOPSIS

   
    
=head1 DESCRIPTION

    This is a module for validating all the phi_base candidate entries and grouping them into a super_hash. 

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <https://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <https://www.ensembl.org/Help/Contact>.

=cut


package Bio::EnsEMBL::EGPipeline::PHIBase::XrefsDbLoader;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::Process');
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::LookUp::LocalLookUp;
use Scalar::Util qw(looks_like_number);
use File::Path qw(make_path remove_tree);

sub param_defaults {

    return {
        'default_branch_code'   => 1,
    };
}

sub fetch_input {
    my $self = shift;
    my $phibase_xrefs = $self->param_required('phibase_xrefs');
    print "------- DB_LOADER ------------     \n";
    
    my $write_is_on = $self->param_required('write');
    if (! $write_is_on) {
      die "Write option is turned off! Changes wont be stored into the DBs.\n
        Run the pipeline with the option -write 1 to overrun this.";
    }

    if (!$phibase_xrefs) {
      die 'phibase_xrefs is not defined'; # Will cause job to fail and leave a message in log_message
    } 

    my $phi_release = $self->param('phi_release');
    if (!$phi_release || $phi_release eq '') {
       die "phi_release parameter missing\n";
    }

    my $division = $self->param('_division');
    my $div = $division;
    my $inputfile  = $self->param('inputfile');


    my ($db_version, $core_host) = $self->_get_params_from_registry();
    $self->param('_db_version',$db_version);
    $self->param('_core_host',$core_host);


    my $output_path = substr $inputfile, 0, rindex( $inputfile, q{/} );
    $output_path = $output_path . "/$div" . "_$db_version" . "_results_PHI_base_$phi_release";
    
    remove_tree($output_path);
    eval { make_path($output_path) };
    if ($@) {
        print "Couldn't create $output_path: $@";
    }
    $self->param('_output_path',$output_path);

    my $out_csv_file =  "$output_path/$div" . "_$db_version" . "_results_$phi_release.csv";
    $self->param('out_file',$out_csv_file);
    
}

=head2 _get_params_from_registry

    Description : Private method that gets the db_version and the mysql_core_host values from the registry file;
    Hoping to find a more elegant way to do this one day...

=cut

sub _get_params_from_registry {
    my $self = shift @_; 

    my $lookup = Bio::EnsEMBL::LookUp::LocalLookUp->new( -SKIP_CONTIGS => 1,
                                                         -NO_CACHE     => 1 );
    
    my $sample_db_name = (keys %{$lookup}->{'dbas'})[0];
    my $version_beg_index = index ($sample_db_name, q{_core})+6;
    my $db_version = substr $sample_db_name, $version_beg_index, 5;
    if (!$db_version) {
        $db_version = 'UNKNOWN_VERSION';
    }
    
    my $sample_tax_id = (keys %{$lookup->{'dbas_by_taxid'}})[0];
    my $dbas = ($lookup->get_all_by_taxon_id ($sample_tax_id));
    my $dba = @$dbas[0];
    my $core_host = $dba->dbc->host();
    if (!$core_host) {
        $db_version = 'UNKNOWN_HOST';
    }

    return ($db_version, $core_host);
}


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_; 

    
    my $phibase_xrefs = $self->param('phibase_xrefs');
    $self->_write_phibase_xrefs($phibase_xrefs);

}


=head2 _write_phibase_xrefs

    Description : Private method that gets values from the phibase_xrefs accumulator and stores them into the DBs;

=cut

sub _write_phibase_xrefs {
    my ($self, $phibase_xrefs) = @_;
    my $phi_release = $self->param('phi_release');
    my $csv_output_file = $self->param('out_file');
    my $filehandler;
    

    if ($csv_output_file) {
        print "trying to open fh\n";
        open($filehandler, '>', $csv_output_file) or die "Could not open file '$csv_output_file' $!";
    } else {
        print "NO FH.  Could not find output file '$csv_output_file'\n";
    }
    
    
    while ( my ( $specific_species_name, $translations ) = each %$phibase_xrefs ) {
        print "............................\nStoring  {$specific_species_name}\n Writing output file to  $csv_output_file\n";

        my $dba = Bio::EnsEMBL::Registry->get_DBAdaptor($specific_species_name, "Core");
        my $dbc = $dba->dbc();
        my $dbname = $dbc->dbname();
        my $dbentry_adaptor = $dba->get_DBEntryAdaptor();
        my $group = 0; # associated_group linking the experiments conditions (phenotype, host and experimental evidence) within a publication/xref
      

        while ( my ( $translation, $phis ) = each %$translations ) {# $translation=translation_id; $phis=nested hash with PHI additional info (phi_number and phenotype, host, experiment info...) 
            while ( my ( $phi, $associated_xrefs ) = each %$phis ) { # $phi= phi id number; $associated_xrefs=rest of PHI info (phenotype, host, experiment refs, conditions and DOI...)
                print "-------------\n{$specific_species_name}{$translation} {$phi} \n";

                my $phi_dbentry =  Bio::EnsEMBL::OntologyXref->new(
                        -PRIMARY_ID  => $phi,
                        -DBNAME      => 'PHI',
                        -DISPLAY_ID  => $phi,
                        -DESCRIPTION => $phi,
                        -RELEASE     => $phi_release,
                        -INFO_TYPE   => 'DIRECT' );

                for my $assoc_xref (@$associated_xrefs) {

                    my $line = $phi . "," . $specific_species_name ;

                    # host species entry
                    my $host_db_entry = Bio::EnsEMBL::DBEntry->new(
                            -PRIMARY_ID => $assoc_xref->{host}{id},
                            -DBNAME     => 'NCBI_TAXONOMY',
                            -RELEASE    => 1,
                            -DISPLAY_ID => $assoc_xref->{host}{label}
                    );
                    $line = $line . "," . $assoc_xref->{host}{id}. ","; 

                    # experimental conditions
                    my @conditions_db_entries; 
                    foreach my $cond ( @{$assoc_xref->{condition}} ) {
            
                        my $condition_db_entry = Bio::EnsEMBL::DBEntry->new(
                            -PRIMARY_ID => $cond->{id},
                            -DBNAME     => 'PHIE',
                            -RELEASE    => 1,
                            -DISPLAY_ID => $cond->{label} 
                        );
                        $line = $line . $cond->{id}. ";"; 
                        push (@conditions_db_entries, $condition_db_entry);
                    }

                    # phenotype entry  
                    my $phenotype_db_entry = Bio::EnsEMBL::DBEntry->new(
                            -PRIMARY_ID => $assoc_xref->{phenotype}{id},
                            -DBNAME     => 'PHIP',
                            -RELEASE    => 1,
                            -DISPLAY_ID => $assoc_xref->{phenotype}{label} # same here, we can spare the repetition with id
                    );
                    $line = $line . "," . $assoc_xref->{phenotype}{id} . "," . $assoc_xref->{evidence}[0] . ","; 

                    # work out literature
                    my $publication_name = 'PUBMED';
                    my $publications     = $assoc_xref->{pubmed};
                    if ( !defined $publications || scalar( @{$publications} ) == 0 ) { # if no PUBMED reference...
                        if ( defined $assoc_xref->{doi} && scalar( @{ $assoc_xref->{doi} } ) > 0 ) { # ...but DOI present
                            $publication_name = 'DOI';
                            $publications     = $assoc_xref->{doi};
                            $line = $line . $assoc_xref->{doi}[0] ;
                        } else {
                            $publications = ['ND'];
                            $line = $line . "ND" ;
                        }
                    } else {
                        foreach my $pubmed ( @{$publications} ) {
                            $line = $line . $pubmed . ";" ;
                        }     
                    }


                    # publications 
                    my $rank = 0; # rank inside the same associated group: 0-phenotype, 1-host, 2-experimental evidence
                    for my $pub (@$publications) {
                        $group++;
                        my $publication_id = lc( _remove_spaces($pub) );
                        my $pub_entry = Bio::EnsEMBL::DBEntry->new(
                            -PRIMARY_ID => $publication_id,
                            -DBNAME     => $publication_name,
                            -DISPLAY_ID => $publication_id,
                            -INFO_TYPE  => $assoc_xref->{evidence}[0] );

                        # add 3 associated_xref per associated_group/publication: phenotype, host, experimental evidence
                        $phi_dbentry->add_associated_xref( 
                            $phenotype_db_entry, 
                            $pub_entry, 
                            'phenotype', 
                            $group, 
                            $rank++ );

                        $phi_dbentry->add_associated_xref( 
                            $host_db_entry,
                            $pub_entry,
                            'host',
                            $group, 
                            $rank++ );  

                        foreach my $condition_db_entry ( @conditions_db_entries ) {
                            $phi_dbentry->add_associated_xref( 
                                $condition_db_entry,
                                $pub_entry, 
                                'experimental evidence',
                                $group, 
                                $rank++ );
                        }
                        $phi_dbentry->add_linkage_type( 'ND', $pub_entry );
                    }
                    
                    my $translation_stable_id = $assoc_xref->{translation_stable_id};
                    my $percent_identity = $assoc_xref->{percent_identity};
                    $line = $line . ",$translation_stable_id,$percent_identity\n";
                    if ($filehandler ) {
                        print $filehandler $line;
                    }
             
                    print $line;
                } ## end for my $assoc_xref (@$associated_xrefs)
                $dbentry_adaptor->store( $phi_dbentry, $translation, 'Translation', 1 );
            }
        }

        print "Applying colors to $specific_species_name\n";
    
        my %gene_color; # phenotype description. 
        $dbc->sql_helper()->execute_no_return(
            -SQL => qq/select t.gene_id, x.display_label from associated_xref a, object_xref o,
             xref x, transcript t, translation tl, seq_region sr, coord_system cs, meta m   
             where a.object_xref_id = o.object_xref_id and condition_type = 'phenotype' 
             and tl.transcript_id=t.transcript_id and x.xref_id = a.xref_id 
             and o.ensembl_id = tl.translation_id and t.seq_region_id=sr.seq_region_id 
             and cs.coord_system_id=sr.coord_system_id and cs.species_id=m.species_id 
             and m.meta_key='species.production_name' and m.meta_value='$specific_species_name';/,
            -CALLBACK => sub {
                my @row   = @{ shift @_ };
                my $color = lc( $row[1] );
                $color =~ s/^\s*(.*)\s*$/$1/;
                $color =~ s/\ /_/g;
                if ( !exists( $gene_color{ $row[0] } ) ) {
                    $gene_color{ $row[0] } = $color;
                } else {
                    if ( $gene_color{ $row[0] } eq $color ) {
                        return;
                    } else {
                        $gene_color{ $row[0] } = 'mixed_outcome';
                    }
                }
                 return;
            } 
        ); 
         # storing PHI-base phenotype of the mutants, and entries for 'PHI' External database in gene_attrib
        foreach my $gene ( keys %gene_color ) {
            print "Setting " . $gene . " as " . $gene_color{$gene} . "\n";
            $dbc->sql_helper()->execute_update(
                -SQL => q/INSERT IGNORE INTO gene_attrib (gene_id, attrib_type_id, value) VALUES ( ?, 358, ?)/,
                -PARAMS => [ $gene, $gene_color{$gene} ] );
            $dbc->sql_helper()->execute_update(
                -SQL => q/INSERT IGNORE INTO gene_attrib (gene_id, attrib_type_id, value) VALUES ( ?, 317, 'PHI')/,
                -PARAMS => [$gene] );
        }   
        $dbc && $dbc->disconnect_if_idle();
    }  
    close $filehandler;
}


sub write_output {
    my ($self, @branch_dbas) = @_;

    my $def_branch_code = $self->param('default_branch_code');
    my $registry_entries = $self->_build_output_hash();
    
    # "fan out" into fdefault_branch_code:
    $self->dataflow_output_id($registry_entries, $def_branch_code);
    
}



=head2 _build_output_hash 

    Description: a private method that returns a hash of parameters for all subjobs. 

=cut

sub _build_output_hash {
    my $self = shift;
    my @hashes = ();
    my $db_version = $self->param_required("_db_version");
    my $output_path = $self->param_required("_output_path");
    my $core_host = $self->param_required("_core_host");
    #                     ----  FIELDS ----
    #                    # '_db_version' # ensembl release version (from registry file)
    #                    # '_core_host' # server with core DBs (from registry file)
    #                    # '_output_path' # path where all output files for this run will be stored   

    my $job_param_hash = {};
    $job_param_hash->{ '_db_version' } = $db_version;
    $job_param_hash->{ '_core_host' } = $core_host;
    $job_param_hash->{ '_output_path' } = $output_path;
    
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
