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

Bio::EnsEMBL::EGPipeline::BlastAlignment::CreateBlastDB

=head1 DESCRIPTION

Convert a fasta file into a BLAST database.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::BlastAlignment::CreateBlastDB;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use File::Basename qw(fileparse);

sub param_defaults {
  return {
    'blast_db'      => undef,
    'database_type' => 'pep',
  };
}

sub fetch_input {
  my ($self) = @_;
  my $makeblastdb_exe = $self->param_required('makeblastdb_exe');
  my $db_fasta_file   = $self->param_required('db_fasta_file');
  
  if (!-e $db_fasta_file) {
    $self->throw("Fasta file '$db_fasta_file' does not exist.");
  }
  
  if (!$self->param_is_defined('blast_db')) {
    my ($name, $dirs, undef) = fileparse($db_fasta_file, qr/\.[^.]*/);
    $self->param('blast_db', "$dirs$name.blastdb");
  }
}

sub run {
  my ($self) = @_;
  my $makeblastdb_exe = $self->param('makeblastdb_exe');
  my $db_fasta_file   = $self->param('db_fasta_file');
  my $blast_db        = $self->param('blast_db');
  my $database_type   = $self->param('database_type');
  
  my $blast_db_type = $database_type eq 'pep' ? 'prot' : 'nucl';
  
  my $cmd = "$makeblastdb_exe -in $db_fasta_file -out $blast_db -dbtype $blast_db_type";
  
  my $out = `$cmd 2>&1`;
  if ($out =~ /error/mi) {
    $self->throw("Error when executing $cmd:\n$out");
  } else {
    my ($suffix_initial) = $blast_db_type =~ /^(\w)/;
    
    my $results = (
      -e "$blast_db.$suffix_initial".'hr' &&
      -e "$blast_db.$suffix_initial".'in' &&
      -e "$blast_db.$suffix_initial".'sq'
    );
    my $results_large_file = (
      -e "$blast_db.00.$suffix_initial".'hr' &&
      -e "$blast_db.00.$suffix_initial".'in' &&
      -e "$blast_db.00.$suffix_initial".'sq' &&
      -e "$blast_db.$suffix_initial".'al'
    );
    
    if (!$results && !$results_large_file) {
      $self->throw("Error when executing $cmd:\n$out");
    }
  }
  
  # To keep the Ensembl Blast runnable happy, need to do this...
  `touch $blast_db` unless -e $blast_db;
}

sub write_output {
  my ($self) = @_;
  my $db_fasta_file     = $self->param_required('db_fasta_file');
  my $proteome_source   = $self->param('proteome_source');
  my $logic_name_prefix = $self->param('logic_name_prefix');
  my $source_species    = $self->param('source_species');
  
  my $output_ids = {
    'blast_db' => $self->param_required('blast_db'),
  };
  
  if (defined $proteome_source && defined $logic_name_prefix) {
    if ($logic_name_prefix eq 'file') {
      my ($name, undef, undef) = fileparse($db_fasta_file, qr/\.[^.]*/);
      ($logic_name_prefix = $name) =~ s/\W/_/g;
    }
    
    if ($proteome_source eq 'database') {
      if ($source_species) {
        $source_species =~ /^(\w)[^_]*_(\w+)/;
        $logic_name_prefix = "$1$2";
      }
    }
    $logic_name_prefix =~ s/\s+/_/g;
    
    $$output_ids{'logic_name_prefix'} = lc($logic_name_prefix);
  }
  
  $self->dataflow_output_id($output_ids, 1);
}

1;
