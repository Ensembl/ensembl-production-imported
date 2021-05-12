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

    Bio::EnsEMBL::EGPipeline::PHIBase::PHISummarySql

=head1 SYNOPSIS

   
=head1 DESCRIPTION

    This is a generic RunnableDB module for producing a summary of all PHI-base associated xrefs that have been stored in the dbs present in the registry for a particular division.
    It's an implementation of ensembl-production's module Bio::EnsEMBL::Production::Pipeline::Common::DbCmd.pm which
    executes a set of SQL commands from a given file by running mysql on the command line.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <https://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <https://www.ensembl.org/Help/Contact>.

=cut


package Bio::EnsEMBL::EGPipeline::PHIBase::PHISummarySql;

use strict;
use warnings;

use File::Path qw(make_path);
use base (
  'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
  'Bio::EnsEMBL::Production::Pipeline::Common::Base',
  'Bio::EnsEMBL::Hive::Process'
);

sub param_defaults {
  my ($self) = @_;
  return {
    %{$self->SUPER::param_defaults},
    'db_type'      => 'core',
    'input_query'  =>    'SELECT /*+ MAX_EXECUTION_TIME(86400000) */ count(DISTINCT g.gene_id) 
                   FROM external_db e join xref x using (external_db_id) 
                   JOIN object_xref ox using (xref_id) 
                   JOIN translation tl on tl.translation_id = ox.ensembl_id 
                   JOIN transcript tc using (transcript_id) 
                   JOIN gene g using(gene_id) 
                   WHERE e.db_name="PHI";


                   SELECT /*+ MAX_EXECUTION_TIME(86400000) */ count(distinct object_xref_id) 
                   FROM external_db e join xref x using (external_db_id) 
                   JOIN object_xref ox using (xref_id) 
                   WHERE e.db_name="PHI" and ox.ensembl_object_type="Translation";

                   SELECT count(distinct species_id)  from external_db e join xref x using (external_db_id)
                   JOIN object_xref ox using (xref_id) 
                   JOIN associated_xref ax using (object_xref_id) 
                   JOIN associated_group ag using (associated_group_id) 
                   JOIN ontology_xref oox using (object_xref_id) 
                   JOIN translation tlt on ox.ensembl_id=tlt.translation_id 
                   JOIN transcript trp using (transcript_id) 
                   JOIN seq_region sr using (seq_region_id) 
                   JOIN coord_system using(coord_system_id) 
                   WHERE e.db_name like "PHI%";',
            
    'append'   => [qw( | grep -v "Warning: Using a password")],
  };
}

sub fetch_input {
  my $self = shift @_;
  my $db_type = $self->param('db_type');
  
  if ($db_type eq 'hive') {
    $self->param('db_conn', $self->dbc);
  } else {
    $self->param('db_conn', $self->get_DBAdaptor($db_type)->dbc);
  }
  
  my $dbname = $self->param_required('dbname');
  my $db_version = $self->param('_db_version');

  my $species = $self->param('species');
  if ($dbname && $dbname ne '' ) {
    print "------------------------------\n SUMMARY of existing phibase xrefs in $dbname \n";
  }

  my $phi_release = $self->param('phi_release');
  if (!$phi_release || $phi_release eq '') {
    die "phi_release parameter missing\n";
  }

  my $general_out_path = $self->param_required("_output_path");
  my $temp_sql_out_path = $general_out_path . "/PHIbase_$db_version" . "_SQL_results";
  eval { make_path($temp_sql_out_path) };
  if ($@) {
    print "Couldn't create $temp_sql_out_path: $@";
  }

  my $temp_file = "$temp_sql_out_path/$dbname.sqr";
  $self->param('output_file',$temp_file);
  print "Writing temp file to $temp_file \n";

  $self->SUPER::fetch_input();
  
}

1;
