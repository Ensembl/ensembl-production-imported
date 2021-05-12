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

Bio::EnsEMBL::EGPipeline::ProjectGeneDesc::EmailReport

=head1 DESCRIPTION

Email a summary of gene descriptions that were projected.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::ProjectGeneDesc::EmailReport;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EmailReport');

use File::Spec::Functions qw(catdir);

sub fetch_input {
  my ($self) = @_;
  
  my $projection_dir = $self->param_required('projection_dir');
  my $store_data     = $self->param_required('store_data');
  my $from_species   = $self->param_required('source');
  my $summary        = $self->param_required('summary');
  
  my $report;
  if (!$store_data) {
    $report .= "Pipeline was run with a flag to prevent storing anything in the database.\n";
  }
  $report .= "Detailed log files are available in: $projection_dir.\n";
  
  my $title = "Projection summary from $from_species";
  my $columns = ['From', 'From genes', 'To', 'To genes', 'Homologies', 'Projected'];
  my $sorted_summary = $self->sort_summary($summary);
  $report .= $self->format_table($title, $columns, $sorted_summary);
  
  my $summary_file = catdir($projection_dir, 'summary.txt');
  open my $data, '>', $summary_file or $self->throw($!);
  print $data $report;
  close $data;
  
  $self->param('text', $report);
}

sub sort_summary {
  my ($self, $summary) = @_;
  
  my %summary;
  foreach my $s (@$summary) {
    my $to = $$s[2];
    $summary{$to} = $s;
  }
  
  my @sorted_summary;
  foreach my $to (sort keys %summary) {
    push @sorted_summary, $summary{$to};
  }
  
  return \@sorted_summary;
}

1;
