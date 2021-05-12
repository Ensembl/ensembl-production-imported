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

Bio::EnsEMBL::EGPipeline::RNAFeatures::EmailRNAGenesReport

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::RNAFeatures::EmailRNAGenesReport;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EmailReport');

sub fetch_input {
  my ($self) = @_;
  my $species             = $self->param_required('species');
  my $use_mirbase         = $self->param_required('use_mirbase');
  my $use_trnascan        = $self->param_required('use_trnascan');
  my $use_cmscan          = $self->param_required('use_cmscan');
  my $mirbase_logic_name  = $self->param_required('mirbase_logic_name');
  my $trnascan_logic_name = $self->param_required('trnascan_logic_name');
  my $cmscan_logic_name   = $self->param_required('cmscan_logic_name');
  
  my $dbh = $self->core_dbh();
  
  my @logic_names;
  push @logic_names, $mirbase_logic_name  if $use_mirbase;
  push @logic_names, $trnascan_logic_name if $use_trnascan;
  push @logic_names, $cmscan_logic_name   if $use_cmscan;
  
  my $reports = "The RNA genes pipeline has completed for $species.\n";
  $reports .= $self->text_summary($dbh, \@logic_names);
  foreach my $logic_name (@logic_names) {
    $reports .= $self->report_one($dbh, $logic_name);
  }
  
  $self->param('text', $reports);
}

sub text_summary {
  my ($self, $dbh, $logic_names) = @_;
  
  my $logic_names_list = join("','", @$logic_names);
  
  my $sql = "
    SELECT 
      logic_name, 
      COUNT(*) AS count_of_genes 
    FROM 
      gene INNER JOIN 
      analysis USING (analysis_id) 
    WHERE logic_name IN ('$logic_names_list') 
    GROUP BY logic_name 
  ;";
  
  my $sth = $dbh->prepare($sql);
  $sth->execute();
  
  my $title = 
    "Genes are linked to the following analyses:";
  my $columns = $sth->{NAME};
  my $results = $sth->fetchall_arrayref();
  
  return $self->format_table($title, $columns, $results);
}

sub report_one {
  my ($self, $dbh, $logic_name) = @_;
  
  my $sql = "
    SELECT 
      biotype, 
      COUNT(*) AS count_of_genes
    FROM  
      gene INNER JOIN 
      analysis USING (analysis_id)
    WHERE logic_name = ? 
    GROUP BY 
      biotype 
    ORDER BY 
      biotype 
  ;";
  
  my $sth = $dbh->prepare($sql);
  $sth->execute($logic_name);
  
  my $title = "Biotype summary ($logic_name):";
  my $columns = $sth->{NAME};
  my $results = $sth->fetchall_arrayref();
  
  return $self->format_table($title, $columns, $results);
}

1;
