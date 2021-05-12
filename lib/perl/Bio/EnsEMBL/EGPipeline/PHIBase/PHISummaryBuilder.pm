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

    Bio::EnsEMBL::EGPipeline::PHIBase::PHISummaryBuilder

=head1 SYNOPSIS
   
    
=head1 DESCRIPTION

    A RunnableDB module for summarizing all results from the sql queries stored by the analysis phi_summary_sql. 
    It reads the different DBs output files with the results from a sql query at a given location and outputs a single file grouping all results together. 

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <https://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <https://www.ensembl.org/Help/Contact>.

=cut


package Bio::EnsEMBL::EGPipeline::PHIBase::PHISummaryBuilder;

use strict;
use warnings;
use Bio::EnsEMBL::LookUp::LocalLookUp;
use Bio::EnsEMBL::Registry;
use Data::Dumper;

use base (
  'Bio::EnsEMBL::Production::Pipeline::Common::Base',
);
use Time::Piece;
use Carp;


sub param_defaults {
  my ($self) = @_;
  
  return {
  };
}

sub fetch_input {
  my $self = shift @_;

  my $db_version = $self->param('_db_version');

  my $general_out_path = $self->param_required("_output_path");
  my $temp_sql_out_path = $general_out_path . "/PHIbase_$db_version" . "_SQL_results";
  $self->param('temp_sql_out_path', $temp_sql_out_path);
  print "------------------------------\n Building phibase xrefs summary from $temp_sql_out_path \n";
}


sub run {
    my $self = shift @_;

    my $temp_sql_out_path = $self->param('temp_sql_out_path'); 
    my $general_out_path = $self->param_required("_output_path");
    my $core_mysql_host = $self->param_required("_core_host");
    my $division = $self->param('_division');
    my $div = $division;# @$division[0];
    my $date = localtime->dmy('_');
    my @data_types = ('GENES','TRANSLATIONS','GENOMES');
    my %global_counts = (
                          'GENES'         => 0,
                          'TRANSLATIONS'  => 0,
                          'GENOMES'       => 0,
                        );

    my $output_file = "$general_out_path/summary_$date.txt";
    my $db_version = $self->param('_db_version');
    open my $output, ">", $output_file or croak "Could not open output_file";
    print "---------- Printing to $output_file\n";

    print $output ">> Total number of EnsemblFungi_$db_version genes and translations in  $core_mysql_host // $date\n";
    my $total_nb_genes = 0;
    my $total_nb_translations = 0;
    
    foreach my $fp (glob("$temp_sql_out_path/*.sqr")) {
      my $db_name = substr $fp, (rindex( $fp, q{/} )+1), rindex ($fp, q{.});
      print $output "---------------------\n$db_name\n";

      open my $fh, "<", $fp or croak "can't read open '$fp'";
      my $dt_ctr = 0;
      while (<$fh>) {
        my $data_type = $data_types[$dt_ctr];
        my $item_count = $_;
        chomp $item_count;
        print $output "$item_count phi annotated $data_type\n" ;
        $global_counts{$data_type} += $item_count;
        $dt_ctr++;
      }
      close $fh or croak "can't read close '$fp'";
    }
    print $output "*******************************\n" ;
    print $output "Total number of Ensembl$div items having PHI-base annotations is " . 
                  $global_counts{'GENES'} .  " genes with " . 
                  $global_counts{'TRANSLATIONS'} . " translations.\n" .
                  "Number of genomes (species or strains) with xrefs in EnsemblFungi: " .
                  $global_counts{'GENOMES'};
    close $output or croak "can't read close '$output'";  
}

1;
