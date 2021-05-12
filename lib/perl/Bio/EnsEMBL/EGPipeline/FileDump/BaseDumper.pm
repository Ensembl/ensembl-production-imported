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

package Bio::EnsEMBL::EGPipeline::FileDump::BaseDumper;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use File::Path qw(make_path);
use File::Spec::Functions qw(catdir);

sub param_defaults {
  return {
    'db_type'            => 'core',
    'gene_centric'       => 0,
    'eg_dir_structure'   => 0,
    'eg_filename_format' => 0,
    'file_varname'       => 'out_file',
  };
}

sub fetch_input {
  my ($self) = @_;
  
  if (defined $self->param('escape_branch') and 
      $self->input_job->retry_count >= $self->input_job->analysis->max_retry_count) 
  {
    $self->dataflow_output_id($self->input_id, $self->param('escape_branch'));
    $self->input_job->autoflow(0);
    $self->complete_early("Failure probably due to memory limit, retrying with a higher limit.");
  }
  
  my $species       = $self->param('species');
  my $results_dir   = $self->param('results_dir');
  my $out_file_stem = $self->param('out_file_stem');
  
  my $out_file;
  if (defined $out_file_stem) {
    $out_file = catdir($results_dir, "$species.$out_file_stem");
  } else {
    $out_file = $self->generate_filename();
  }
  
  $self->param('out_file', $out_file);
  $self->param('out_files', [$out_file]);
}

sub write_output {
  my ($self) = @_;
  my $file_varname = $self->param_required('file_varname');
  
  foreach my $out_file (@{$self->param('out_files')}) {
    $self->dataflow_output_id({$file_varname => $out_file}, 1);
  }
}

sub generate_filename {
  my ($self) = @_;
  
  my $species            = $self->param('species');
  my $file_type          = $self->param('file_type');
  my $results_dir        = $self->param('results_dir');
  my $filename           = $self->param('filename');
  my $eg_dir_structure   = $self->param('eg_dir_structure');
  my $eg_filename_format = $self->param('eg_filename_format');
  
  if ($eg_dir_structure) {
    my ($division, $collection) = $self->get_division();
    $results_dir = catdir($results_dir, $division, $file_type, $collection, $species);
    $self->param('results_dir', $results_dir);
  }

  make_path($results_dir);
  
  if (!$filename) {
    if ($eg_filename_format) {
      $filename = $self->generate_eg_filename($species, $file_type);
    } else {
      $filename = $self->generate_vb_filename($species, $file_type);
    }
  }
  
  return catdir($results_dir, $filename);
}

sub generate_eg_filename {
  my ($self, $species, $file_type) = @_;
  
  my $dba = $self->core_dba;
  my $dbname = $dba->dbc->dbname();
  my $assembly = $dba->get_MetaContainer()->single_value_by_key('assembly.default');
  my ($eg_version) = $dbname =~ /([^_]+)_[^_]+_[^_]+$/;
  my $filename = ucfirst($species).".$assembly.$eg_version.$file_type";
  
  return $filename;
}

sub generate_vb_filename {
  my ($self, $species, $file_type) = @_;
  
  my $data_type    = $self->param('data_type');
  my $gene_centric = $self->param('gene_centric');
  
  my $dba = $self->core_dba;
  
  $species =~ s/_/-/;
  $species =~ s/_.+$//;
  
  my $strain  = $dba->get_MetaContainer()->single_value_by_key('species.strain');
  $strain =~ s/[\s\/]+/\-/g;
  
  my $version;
  if ($gene_centric) {
    $version = $dba->get_MetaContainer()->single_value_by_key('genebuild.version');
  } else {
    $version = $dba->get_MetaContainer()->single_value_by_key('assembly.default');
  }
  
  my $filename = ucfirst($species).'-'.$strain.'_'.uc($data_type).'_'."$version.$file_type";
  
  $dba->dbc && $dba->dbc->disconnect_if_idle();
  
  return $filename;
}

sub get_division {
  my ($self) = @_;
  
  my $dba = $self->core_dba;  
  my ($division, $collection);
  if ($dba->dbc->dbname() =~ /(\w+)(\_[0-9a-zA-Z]+_collection)_/) {
    $division = $1;
    $collection = "$division$2";
  } else {
    $division = $dba->get_MetaContainer->get_division();
    $division = lc($division);
    $division =~ s/ensembl//;
  }
  return ($division, $collection);
}

1;
