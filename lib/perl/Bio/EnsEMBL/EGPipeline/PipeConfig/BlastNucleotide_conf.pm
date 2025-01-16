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

Bio::EnsEMBL::EGPipeline::PipeConfig::BlastNucleotide_conf

=head1 DESCRIPTION

Align nucleotide sequences against a set of core databases.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::PipeConfig::BlastNucleotide_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf');

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use Bio::EnsEMBL::Hive::Version v2.4;

use File::Spec::Functions qw(catdir);

sub default_options {
  my $self = shift @_;
  return {
    %{$self->SUPER::default_options},
    
    pipeline_name => 'blast_nucleotide_'.$self->o('ensembl_release'),
    
    species      => [],
    antispecies  => [],
    division     => [],
    run_all      => 0,
    meta_filters => {},
    
    # Parameters for dumping and splitting Fasta query files.
    max_seq_length          => 10000000,
    max_seq_length_per_file => $self->o('max_seq_length'),
    max_seqs_per_file       => 500,
    max_files_per_directory => 100,
    max_dirs_per_directory  => $self->o('max_files_per_directory'),
    
    max_hive_capacity => 100,
    
    makeblastdb_exe  => 'makeblastdb',
    blast_exe        => 'blastn',
    blast_threads    => 3,
    blast_parameters => '-word_size 11 -num_alignments 100000 -num_descriptions 100000 -lcase_masking -num_threads '.$self->o('blast_threads'),
    
    # For parsing the output.
    output_regex     => '^\s*(\w+)',
    pvalue_threshold => 0.01,
    filter_prune     => 1,
    filter_min_score => 200,
    
    # Specify blast_db if you don't want it in the directory alongside db_file.
    blast_db => undef,
    
    # Filter so that only best X hits get saved.
    blast_top_x => 1,
    
    # blast results go in a core db by default.
    db_type => 'core',
    
    # Remove existing *_align_features; if => 0 then existing analyses
    # and their features will remain, with the logic_name suffixed by '_bkp'.
    delete_existing => 1,

    # Retrieve analysis descriptions from the production database;
    # the supplied registry file will need the relevant server details.
    production_lookup => 1,

  };
}

sub beekeeper_extra_cmdline_options {
  my ($self) = @_;
  
  my $options = join(' ',
    $self->SUPER::beekeeper_extra_cmdline_options,
    "-reg_conf ".$self->o('registry')
  );
  
  return $options;
}

sub hive_meta_table {
  my ($self) = @_;
  
  return {
    %{$self->SUPER::hive_meta_table},
    'hive_use_param_stack'  => 1,
  };
}

sub pipeline_create_commands {
  my ($self) = @_;
  
  return [
    @{$self->SUPER::pipeline_create_commands},
    'mkdir -p '.$self->o('pipeline_dir'),
  ];
}

sub pipeline_wide_parameters {
 my ($self) = @_;
 
 return {
   %{$self->SUPER::pipeline_wide_parameters},
   'db_type' => $self->o('db_type'),
 };
}

sub pipeline_analyses {
  my $self = shift @_;
  
  return [
    {
      -logic_name      => 'CreateBlastDB',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::CreateBlastDB',
      -max_retry_count => 0,
      -input_ids       => [{}],
      -parameters      => {
                            makeblastdb_exe   => $self->o('makeblastdb_exe'),
                            db_fasta_file     => $self->o('db_fasta_file'),
                            blast_db          => $self->o('blast_db'),
                            database_type     => 'nucl',
                          },
      -rc_name         => 'normal',
      -flow_into       => ['SpeciesFactory'],
    },

    {
      -logic_name      => 'SpeciesFactory',
      -module          => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EGSpeciesFactory',
      -max_retry_count => 1,
      -parameters      => {
                            species         => $self->o('species'),
                            antispecies     => $self->o('antispecies'),
                            division        => $self->o('division'),
                            run_all         => $self->o('run_all'),
                            meta_filters    => $self->o('meta_filters'),
                            chromosome_flow => 0,
                            regulation_flow => 0,
                            variation_flow  => 0,
                          },
      -rc_name         => 'normal',
      -flow_into       => {
                            '2->A' => WHEN(
                                        '#db_type# eq "core"'          => ['BackupCoreDatabase'],
                                        '#db_type# eq "otherfeatures"' => ['CheckOFDatabase'],
                                      ),
                            'A->2' => ['DumpGenome'],
                          },
      -meadow_type     => 'LOCAL',
    },

    {
      -logic_name      => 'BackupCoreDatabase',
      -module          => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::DatabaseDumper',
      -can_be_empty      => 1,
      -max_retry_count => 1,
      -parameters      => {
                            output_file => catdir($self->o('pipeline_dir'), '#species#', 'core_bkp.sql.gz'),
                          },
      -rc_name         => 'normal',
    },

    {
      -logic_name      => 'CheckOFDatabase',
      -module          => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::CheckOFDatabase',
      -can_be_empty    => 1,
      -max_retry_count => 1,
      -parameters      => {},
      -rc_name         => 'normal',
      -flow_into       => {
                            '2' => ['BackupOFDatabase'],
                            '3' => ['CreateOFDatabase'],
                          },
      -meadow_type     => 'LOCAL',
    },

    {
      -logic_name      => 'BackupOFDatabase',
      -module          => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::DatabaseDumper',
      -can_be_empty    => 1,
      -max_retry_count => 1,
      -parameters      => {
                            db_type     => 'otherfeatures',
                            output_file => catdir($self->o('pipeline_dir'), '#species#', 'otherfeatures_bkp.sql.gz'),
                          },
      -rc_name         => 'normal',
    },

    {
      -logic_name      => 'CreateOFDatabase',
      -module          => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::CreateOFDatabase',
      -can_be_empty    => 1,
      -max_retry_count => 1,
      -parameters      => {},
      -rc_name         => 'normal',
      -flow_into       => ['BackupOFDatabase'],
    },

    {
      -logic_name        => 'DumpGenome',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::DumpGenome',
      -can_be_empty      => 1,
      -analysis_capacity => 5,
      -parameters        => {
                              genome_dir => catdir($self->o('pipeline_dir'), '#species#', 'genome'),
                            },
      -rc_name           => 'normal',
      -flow_into         => {
                              '1->A' => ['AnalysisSetup'],
                              'A->1' => ['SplitGenome'],
                            },
    },

    {
      -logic_name      => 'AnalysisSetup',
      -module          => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::AnalysisSetup',
      -can_be_empty    => 1,
      -max_retry_count => 0,
      -batch_size      => 10,
      -parameters      => {
                            logic_name         => $self->o('logic_name'),
                            display_label      => 'Nucleotide alignments',
                            description        => 'Nucleotide sequences aligned to the genome with <em>blastn</em>.',
                            displayable        => 1,
                            web_data           => '{"type" => "cdna"}',
                            db                 => $self->o('db_fasta_file'),
                            program            => 'blastn',
                            program_file       => $self->o('blast_exe'),
                            parameters         => $self->o('blast_parameters'),
                            module             => 'Bio::EnsEMBL::Analysis::Runnable::BlastEG',
                            linked_tables      => ['dna_align_feature'],
                            db_type            => $self->o('db_type'),
                            db_backup_required => 1,
                            db_backup_file     => catdir($self->o('pipeline_dir'), '#species#', '#db_type#_bkp.sql.gz'),
                            delete_existing    => $self->o('delete_existing'),
                            production_lookup  => $self->o('production_lookup'),
                            production_db      => $self->o('production_db'),
                          },
      -meadow_type     => 'LOCAL',
    },

    {
      -logic_name        => 'SplitGenome',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::FastaSplit',
      -can_be_empty      => 1,
      -analysis_capacity => 5,
      -parameters        => {
                              fasta_file              => '#genome_file#',
                              max_seq_length_per_file => $self->o('max_seq_length_per_file'),
                              max_seqs_per_file       => $self->o('max_seqs_per_file'),
                              max_files_per_directory => $self->o('max_files_per_directory'),
                              max_dirs_per_directory  => $self->o('max_dirs_per_directory'),
                            },
      -rc_name           => 'normal',
      -flow_into         => {
                              '2->A' => ['BlastFactory'],
                              'A->1' => ['FilterHits'],
                            },
    },

    {
      -logic_name      => 'BlastFactory',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::BlastFactory',
      -can_be_empty    => 1,
      -hive_capacity   => $self->o('max_hive_capacity'),
      -max_retry_count => 1,
      -parameters      => {
                            max_seq_length => $self->o('max_seq_length'),
                            queryfile      => '#split_file#',
                          },
      -rc_name         => 'normal',
      -flow_into       => {
                            '2' => ['Blast'],
                          },
    },

    {
      -logic_name      => 'Blast',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::Blast',
      -can_be_empty    => 1,
      -hive_capacity   => $self->o('max_hive_capacity'),
      -max_retry_count => 1,
      -parameters      => {
                            logic_name       => $self->o('logic_name'),
                            query_type       => 'dna',
                            database_type    => 'dna',
                            output_regex     => $self->o('output_regex'),
                            pvalue_threshold => $self->o('pvalue_threshold'),
                            filter_prune     => $self->o('filter_prune'),
                            filter_min_score => $self->o('filter_min_score'),
                            escape_branch    => -1,
                          },
      -rc_name         => '8GB_threads',
      -flow_into       => {
                            '-1' => ['Blast_HighMem'],
                          },
    },

    {
      -logic_name      => 'Blast_HighMem',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::Blast',
      -can_be_empty    => 1,
      -hive_capacity   => $self->o('max_hive_capacity'),
      -max_retry_count => 1,
      -parameters      => {
                            logic_name       => $self->o('logic_name'),
                            query_type       => 'dna',
                            database_type    => 'dna',
                            output_regex     => $self->o('output_regex'),
                            pvalue_threshold => $self->o('pvalue_threshold'),
                            filter_prune     => $self->o('filter_prune'),
                            filter_min_score => $self->o('filter_min_score'),
                            escape_branch    => -1,
                          },
      -rc_name         => '16GB_threads',
      -flow_into       => {
                            '-1' => ['Blast_HigherMem'],
                          },
    },

    {
      -logic_name      => 'Blast_HigherMem',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::Blast',
      -can_be_empty    => 1,
      -hive_capacity   => $self->o('max_hive_capacity'),
      -max_retry_count => 1,
      -parameters      => {
                            logic_name       => $self->o('logic_name'),
                            query_type       => 'dna',
                            database_type    => 'dna',
                            output_regex     => $self->o('output_regex'),
                            pvalue_threshold => $self->o('pvalue_threshold'),
                            filter_prune     => $self->o('filter_prune'),
                            filter_min_score => $self->o('filter_min_score'),
                          },
      -rc_name         => '32GB_threads',
    },

    {
      -logic_name      => 'FilterHits',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::FilterHits',
      -can_be_empty    => 1,
      -max_retry_count => 1,
      -parameters      => {
                            filter_top_x => $self->o('blast_top_x'),
                            logic_name   => $self->o('logic_name'),
                          },
      -rc_name         => 'normal',
    },

  ];
}

sub resource_classes {
  my ($self) = @_;
  
  my $blast_threads = $self->o('blast_threads');

  return {
    %{$self->SUPER::resource_classes},
    '8GB_threads' => {'LSF' => '-q ' . $self->o('queue_name') . ' -n ' . ($blast_threads + 1) . ' -M 8000 -R "rusage[mem=8000,tmp=8000]"'},
    '16GB_threads' => {'LSF' => '-q ' . $self->o('queue_name') . ' -n ' . ($blast_threads + 1) . ' -M 16000 -R "rusage[mem=16000,tmp=16000]"'},
    '32GB_threads' => {'LSF' => '-q ' . $self->o('queue_name') . ' -n ' . ($blast_threads + 1) . ' -M 32000 -R "rusage[mem=32000,tmp=32000]"'},
  }
}

1;
