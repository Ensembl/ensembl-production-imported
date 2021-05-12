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

package Bio::EnsEMBL::EGPipeline::BlastAlignment::AnalysisFactory;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub param_defaults {
  return {};
}

sub run {
  my ($self) = @_;
  my $analyses          = $self->param_required('analyses');
  my $logic_name_prefix = $self->param_required('logic_name_prefix');
  my $db_file           = $self->param_required('db_fasta_file');
  my $blastp            = $self->param_required('blastp');
  my $blastx            = $self->param_required('blastx');
  my $source_species    = $self->param('source_species');
  
  my $proteomic_analyses = [];
  my $genomic_analyses = [];
  
  foreach my $analysis (@{$analyses}) {
    $$analysis{'logic_name'} = $logic_name_prefix.'_'.$$analysis{'logic_name'};
    $$analysis{'db_file'} = $db_file;
    
    if ($source_species) {
      $self->param('species', $source_species);
      my $mc = $self->core_dba->get_MetaContainer();
      my $db_version   = $mc->single_value_by_key('genebuild.version');
      my $species_name = $mc->single_value_by_key('species.scientific_name');
      
      $$analysis{'db_version'} = $db_version;
      $$analysis{'description'} .= " Source database: <em>$species_name</em> proteome (version $db_version).";
      $$analysis{'display_label'} = "$db_version proteins";
    } else {
      $db_file =~ s!.*/([^/]+)$!$1!;
      $$analysis{'description'} .= " Source database: $db_file.";
    }
    
    my $top_x;
    if ($$analysis{'logic_name'} =~ /blastp$/) {
      $top_x = $self->param('blastp_top_x');
    } elsif ($$analysis{'logic_name'} =~ /blastx$/) {
      $top_x = $self->param('blastx_top_x');
    }
    if (defined $top_x) {
      if ($top_x eq 1) {
        $$analysis{'description'} .= " Results are filtered to show the best unique hit, based on E-value.";
      } else {
        $$analysis{'description'} .= " Results are filtered to show the top $top_x unique hits, based on E-value.";
        $$analysis{'description'} .= " Multiple hits in the same region are automatically joined into a single feature.";
      }
    }
    
    if ($$analysis{'logic_name'} =~ /blastp$/ && $blastp) {
      push @$proteomic_analyses, $analysis;
    } elsif ($$analysis{'logic_name'} =~ /blastx$/ && $blastx) {
      push @$genomic_analyses, $analysis;
    }
  }
  
  $self->param('proteomic_analyses', $proteomic_analyses);
  $self->param('genomic_analyses', $genomic_analyses);
}

sub write_output {
  my ($self) = @_;
  my $blastp = $self->param_required('blastp');
  my $blastx = $self->param_required('blastx');
  
  if ($blastp) {
    $self->dataflow_output_id($self->param('proteomic_analyses'), 2);
    $self->dataflow_output_id({}, 4);
  }
  
  if ($blastx) {  
    $self->dataflow_output_id($self->param('genomic_analyses'), 3);
    $self->dataflow_output_id({}, 5);
  }
}

1;
