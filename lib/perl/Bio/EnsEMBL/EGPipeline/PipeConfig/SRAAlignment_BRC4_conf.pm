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

package Bio::EnsEMBL::EGPipeline::PipeConfig::SRAAlignment_BRC4_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf');

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use Bio::EnsEMBL::Hive::Version 2.5;

use File::Spec::Functions qw(catdir);

sub default_options {
  my ($self) = @_;
  return {
    %{ $self->SUPER::default_options() },

    # Pipeline parameters
    pipeline_name => 'sra_alignment_brc4_'.$self->o('ensembl_release'),
    pipeline_dir => $self->o('pipeline_dir'),
    index_dir => catdir($self->o('pipeline_dir'), 'index'),

    # INPUT
    # Datasets to align
    datasets_file => $self->o('datasets_file'),

    # Species factory
    species => [],
    antispecies => [],
    division => [],
    run_all => 0,
    meta_filters => {},
    
    # BEHAVIOUR
    # Use DNA-Seq specific parameters
    # In that case, do not remove bam files
    dnaseq => 0,
    skip_cleanup => $self->o('dnaseq'),
    
    # Download via SRA from NCBI if ENA fails (fallback)
    use_ncbi => 0,
    sra_dir => undef,
    
    # Clean up temp files at the end (split bam etc.)
    clean_up      => 1,
    
    # If there is already an alignment for a sample, redo its htseq-count only
    # (do nothing otherwise)
    redo_htseqcount => 0,
    features => ['exon', 'CDS'],
    alignments_dir => undef,
    
    # Use input metadata, instead of inferring them (pair/strand)
    infer_metadata => 1,

    # For heavy analyses, use multiple cpus
    threads    => 4,
    samtobam_memory => 16000,

    ###########################################################################
    # PATHS
    # Path to trimmomatic binary and adapters folders
    trimmomatic_bin => undef,
    trim_adapters_pe => undef,
    trim_adapters_se => undef,
    
    # Path to hisat2 binaries
    hisat2_dir    => undef,
    
    # If a directory is not specified then the version in the
    # Ensembl software environment will be used.
    samtools_dir  => undef,
    bedtools_dir  => undef,
    bcftools_dir  => undef,
    ucscutils_dir => undef,
    bamutils_dir  => undef,

    ###########################################################################
    # PARAMETERS unlikely to be changed
    
    # Dump genome: repeatmasking
    repeat_masking     => 'soft',
    repeat_logic_names => [],
    min_slice_length   => 0,

    # Aligner options
    max_intron => 1,
  };
}

# Force an automatic loading of the registry in all workers.
sub beekeeper_extra_cmdline_options {
  my ($self) = @_;

  my $options = join(' ',
    $self->SUPER::beekeeper_extra_cmdline_options,
    "-reg_conf ".$self->o('registry')
  );
  
  return $options;
}

sub pipeline_create_commands {
  my ($self) = @_;
  
  # To store alignment commands
  my $align_cmds_table =
    'CREATE TABLE align_cmds ('.
      'auto_id INT AUTO_INCREMENT PRIMARY KEY, '.
      'sample_name varchar(255) NULL, '.
      'cmds text NOT NULL, '.
      'version varchar(255) NULL)';

  return [
    @{$self->SUPER::pipeline_create_commands},
    'mkdir -p '.$self->o('pipeline_dir'),
    'mkdir -p '.$self->o('results_dir'),
    $self->db_cmd($align_cmds_table),
  ];
}

sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{$self->SUPER::pipeline_wide_parameters},
    'redo_htseqcount'      => $self->o('redo_htseqcount'),
    'infer_metadata' => $self->o('infer_metadata'),
    'dnaseq'      => $self->o('dnaseq'),
  };
}

sub pipeline_analyses {
  my ($self) = @_;

  my $aligner = 'hisat2';
  my $aligner_class = 'Bio::EnsEMBL::EGPipeline::Common::Aligner::HISAT2Aligner';
  my $aligner_dir = $self->o('hisat2_dir');

  my @backbone = (
    {
      -logic_name => 'Start',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -input_ids  => [{}],
      -flow_into  => {
        '1->A' => 'Species_factory',
        'A->1' => 'Email_report',
      },
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
      -analysis_capacity => 1,
    },

    {
      -logic_name => 'Species_factory',
      -module     => 'Bio::EnsEMBL::Production::Pipeline::Common::SpeciesFactory',
      -parameters => {
        run_all      => 1,
      },
      -flow_into  => { 2 => 'Species_filter' },
      -max_retry_count => 0,
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
      -analysis_capacity => 1,
    },

    {
      -logic_name => 'Species_filter',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::DatasetSpeciesFactory',
      -parameters => {
        datasets_file => $self->o('datasets_file'),
      },
      -flow_into  => {
        '2->A' => 'Prepare_genome',
        'A->2' => 'Species_report',
      },
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
      -analysis_capacity => 1,
      -max_retry_count => 0,
    },

    {
      -logic_name        => 'Species_report',
      #-module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::EmailReport',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -parameters        => {
                              email        => $self->o('email'),
                              subject      => 'Short Read Alignment pipeline: Report for #species#',
                              samtools_dir => $self->o('samtools_dir'),
                            },
      -rc_name           => 'normal',
      -analysis_capacity => 1,
      -max_retry_count => 0,
    },

    {
      -logic_name        => 'Email_report',
      #-module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::EmailReport',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -parameters        => {
                              email        => $self->o('email'),
                              subject      => 'Short Read Alignment pipeline: final report',
                            },
      -rc_name           => 'normal',
      -analysis_capacity => 1,
      -max_retry_count => 0,
    },

    # Prepare genome
    {
      -logic_name        => 'Prepare_genome',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::PrepareGenome',
      -parameters        => {
        species => '#species#',
        index_dir => $self->o('index_dir'),
        results_dir => $self->o('results_dir'),
        pipeline_dir => $self->o('pipeline_dir'),
      },
      -analysis_capacity => 1,
      -max_retry_count => 0,
      -rc_name           => 'normal',
      -flow_into         => {
        1 => { 'Genome_data_factory' => INPUT_PLUS() },
      },
    },

    {
      -logic_name        => 'Genome_data_factory',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name           => 'normal',
      -flow_into         => {
        '1->A' => { 'Genome_factory' => INPUT_PLUS() },
        'A->1' => { 'Runs_factory' => INPUT_PLUS() },
      },
    },
  );

    ####################################################################
    # Genomes

  my @genomes = (
    {
      -logic_name        => 'Genome_factory',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'mkdir -p #genome_dir#',
      },
      -rc_name           => 'normal',
      -flow_into         => 'Dump',
      -max_retry_count => 0,
    },

    {
      -logic_name        => 'Dump',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name           => 'normal',
      -flow_into         => [
       'GenesCheck',
        WHEN('not -s #genome_file#', 'Dump_genome'),
      ],
    },

    {
      -logic_name        => 'GenesCheck',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::GenesCheck',
      -analysis_capacity => 1,
      -max_retry_count => 0,
      -rc_name           => 'normal',
      -flow_into         => {
        2 => [
          WHEN('not -s #genome_bed_file#', 'Dump_gff3'),
          WHEN('not -s #genome_gtf_file#', 'Dump_gtf'),
        ]
      },
    },

    {
      -logic_name        => 'Dump_gff3',
      -module            => 'Bio::EnsEMBL::EGPipeline::FileDump::GFF3Dumper',
      -analysis_capacity => 10,
      -can_be_empty      => 1,
      -max_retry_count   => 1,
      -parameters        => {
        data_type          => 'basefeatures',
        feature_type       => ['Exon'],
        results_dir        => '#genome_dir#',
      },
      -rc_name           => 'normal',
      -flow_into         => {
        '1' => 'Rename_gff',
      },
    },
    {
      -logic_name        => 'Rename_gff',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'mv "#out_file#" #genome_gff_file#',
      },
      -rc_name           => 'normal',
      -flow_into         => 'gff3_to_bed',
      -max_retry_count => 0,
    },

    {
      -logic_name        => 'gff3_to_bed',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::GFF3_to_bed',
      -rc_name           => 'normal',
      -max_retry_count => 0,
      -parameters        => {
        gff_file          => '#genome_gff_file#',
        bed_file          => '#genome_bed_file#',
      },
    },

    {
      -logic_name        => 'Dump_gtf',
      -module            => 'Bio::EnsEMBL::EGPipeline::FileDump::GTFDumper',
      -analysis_capacity => 10,
      -can_be_empty      => 1,
      -max_retry_count   => 0,
      -parameters        => {
        results_dir => '#genome_dir#',
      },
      -rc_name           => 'normal',
      -flow_into => 'Rename_gtf',
    },
    {
      -logic_name        => 'Rename_gtf',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'mv "#out_file#" #genome_gtf_file#',
      },
      -rc_name           => 'normal',
      -max_retry_count => 0,
    },

    {
      -logic_name        => 'Dump_genome',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::DumpGenome',
      -analysis_capacity => 5,
      -max_retry_count => 0,
      -batch_size        => 2,
      -parameters        => {
        genome_dir         => '#genome_dir#',
        repeat_masking     => $self->o('repeat_masking'),
        repeat_logic_names => $self->o('repeat_logic_names'),
        min_slice_length   => $self->o('min_slice_length'),
      },
      -rc_name           => 'normal',
      -flow_into         => ['Index_genome', 'Sequences_lengths'],
    },

    {
      -logic_name        => 'Index_genome',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::IndexGenome',
      -analysis_capacity => 5,
      -max_retry_count => 0,
      -batch_size        => 2,
      -parameters        => {
        aligner_class => $aligner_class,
        aligner_dir   => $aligner_dir,
        samtools_dir  => $self->o('samtools_dir'),
        memory_mode   => 'default',
        overwrite     => 0,
        escape_branch => -1,
        threads       => $self->o('threads'),
      },
      -rc_name           => '8GB_multicpu',
      -flow_into         => {
        '-1' => 'Index_genome_HighMem',
      },
    },

    {
      -logic_name        => 'Index_genome_HighMem',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::IndexGenome',
      -analysis_capacity => 5,
      -max_retry_count => 0,
      -batch_size        => 2,
      -can_be_empty      => 1,
      -parameters        => {
        aligner_class => $aligner_class,
        aligner_dir   => $aligner_dir,
        samtools_dir  => $self->o('samtools_dir'),
        memory_mode   => 'himem',
        threads       => $self->o('threads'),
      },
      -rc_name           => '32GB_multicpu',
    },

    {
      -logic_name        => 'Sequences_lengths',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::SequenceLengths',
      -analysis_capacity => 5,
      -max_retry_count => 0,
      -batch_size        => 2,
      -can_be_empty      => 1,
      -parameters        => {
        fasta_file  => '#genome_file#',
        length_file => '#length_file#',
      },
      -rc_name           => 'normal',
    },
  );

    ####################################################################
    # Groups and files

  my @files = (
    # Groups
    {
      -logic_name        => 'Runs_factory',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::RunFactory',
      -parameters        => {
        results_dir => '#species_results_dir#',
        datasets_file => $self->o('datasets_file'),
        tax_id_restrict => 0,
        check_library => 0,
      },
      -analysis_capacity => 1,
      -max_retry_count => 0,
      -rc_name           => 'normal',
      -flow_into         => {
        '2->A' => WHEN("not #redo_htseqcount#", 'SRASeqFileFromENA'),
        'A->4' => WHEN("not #redo_htseqcount#", 'Merge_fastq'),
        '3' => WHEN("#redo_htseqcount#", 'GenesCheckForRedo'),
      },
    },

    # HTseq recount path
    {
      -logic_name        => 'GenesCheckForRedo',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::GenesCheck',
      -analysis_capacity => 1,
      -max_retry_count => 0,
      -rc_name           => 'normal',
      -flow_into         => {
        2 => { "SyncAlignmentFiles" => INPUT_PLUS() },
      },
    },

    {
      -logic_name        => 'SyncAlignmentFiles',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::SyncAlignmentFiles',
      -analysis_capacity => 4,
      -max_retry_count   => 0,
      -rc_name           => 'datamove',
      -parameters => {
        alignments_dir => $self->o('alignments_dir'),
      },
      -flow_into         => {
        2 => "ReadMetadata",
      },
    },

    {
      -logic_name        => 'ReadMetadata',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::ReadMetadata',
      -parameters        => {
                              results_dir => '#sample_dir#',
                            },
      -rc_name           => 'normal',
      -analysis_capacity => 1,
      -max_retry_count => 0,
      -flow_into         => {
        2 => { 'HtseqFactory_redo' => INPUT_PLUS() },
      },

    },
    {
      -logic_name        => 'HtseqFactory_redo',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::HtseqFactory',
      -parameters        => {
        bam_file => '#sorted_bam_file#',
        features => $self->o('features'),
      },
      -rc_name           => 'normal',
      -analysis_capacity => 1,
      -max_retry_count => 0,
      -flow_into         => {
        2 => { 'HtseqCount_redo' => INPUT_PLUS() },
      },
    },

    {
      -logic_name        => 'HtseqCount_redo',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::HtseqCountName',
      -parameters        => {
        gtf_file       => '#genome_gtf_file#',
        results_dir    => '#sample_dir#',
      },
      -rc_name           => 'normal',
      -analysis_capacity => 25,
      -max_retry_count => 0,
    },

    # Normal alignment path
    {
      -logic_name        => 'Merge_fastq',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::MergeFastq',
      -parameters        => {
        work_dir => '#species_work_dir#',
      },
      -analysis_capacity => 1,
      -max_retry_count => 0,
      -rc_name           => 'normal',
      -flow_into         => {
        2 => { "PreAlignment" => INPUT_PLUS() },
      }
    },

    {
      -logic_name        => 'SRASeqFileFromENA',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::SRASeqFile',
      -analysis_capacity => 4,
      -max_retry_count => 3,
      -failed_job_tolerance => 10,
      -parameters        => {
        work_dir => '#species_work_dir#',
        use_ncbi => $self->o('use_ncbi'),
      },
      -rc_name           => 'normal',
      -flow_into         => {
        2 => 'ConvertSpace',
        3 => 'SRASeqFileFromNCBI',
      },
    },

    {
      -logic_name        => 'SRASeqFileFromNCBI',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::SRASeqFileFromNCBI',
      -analysis_capacity => 4,
      -max_retry_count => 2,
      -can_be_empty      => 1,
      -parameters        => {
        work_dir => '#species_work_dir#',
        sra_dir => $self->o('sra_dir'),
      },
      -flow_into         => { 2 => 'ConvertSpace' },
      -rc_name           => 'normal',
    },

    {
      -logic_name        => 'ConvertSpace',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::ConvertSpace',
      -analysis_capacity => 4,
      -max_retry_count => 0,
      -parameters        => {
        seq_file_1 => '#run_seq_file_1#',
        seq_file_2 => '#run_seq_file_2#',
      },
      -flow_into         => {
        2 => [
          '?accu_name=run_seq_files_1&accu_address=[]&accu_input_variable=run_seq_file_1',
          '?accu_name=run_seq_files_2&accu_address=[]&accu_input_variable=run_seq_file_2',
        ],
      },
      -rc_name           => 'normal',
    },

    {
      -logic_name        => 'PreAlignment',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name           => 'normal',
      -analysis_capacity => 1,
      -flow_into         => {
        '1' => WHEN('#trim_reads#', { 'Trim' => INPUT_PLUS() }, ELSE { 'CheckMetadata' => INPUT_PLUS() })
      },
    },

    {
      -logic_name        => 'Trim',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::Trim',
      -analysis_capacity => 4,
      -max_retry_count => 0,
      -can_be_empty      => 1,
      -parameters        => {
        seq_file_1     => '#sample_seq_file_1#',
        seq_file_2     => '#sample_seq_file_2#',
        threads        => $self->o('threads'),
        trimmomatic_bin => $self->o('trimmomatic_bin'),
        trim_adapters_pe => $self->o('trim_adapters_pe'),
        trim_adapters_se => $self->o('trim_adapters_se'),
      },
      -rc_name           => '16GB',
      -flow_into         => {
        '2' => 'CheckMetadata',
      },
    },
  );

    ####################################################################
    # Alignment
  my @sample_alignment = (

    {
      -logic_name        => 'CheckMetadata',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::GenesCheck',
      -analysis_capacity => 1,
      -max_retry_count => 0,
      -rc_name           => 'normal',
      -flow_into         => {
        2 => WHEN("#infer_metadata#", { "SubsetSequence" => INPUT_PLUS() }, ELSE({ "UseInputMetadata" => INPUT_PLUS() })),
        3 => { "UseInputMetadata" => INPUT_PLUS() },
      },
    },

    {
      -logic_name        => 'UseInputMetadata',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name           => 'normal',
      -analysis_capacity => 1,
      -max_retry_count => 0,
      -flow_into         => {
        '1' => { "AlignSequence" => {
            is_paired => "#input_is_paired#",
            is_stranded => "#input_is_stranded#",
            strand_direction => "#input_strand_direction#",
            strandness => "#input_strandness#",
          } }
      },
    },
    
    # Subset alignment to infer strand-specificity
    {
      -logic_name        => 'SubsetSequence',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::FastqSubset',
      -parameters        => {
        seq_file_1     => '#sample_seq_file_1#',
        seq_file_2     => '#sample_seq_file_2#',
      },
      -analysis_capacity => 50,
      -max_retry_count => 0,
      -rc_name           => 'normal',
      -flow_into         => {
        '1' => { 'AlignSubsetSequence' => INPUT_PLUS() },
      },
    },

    {
      -logic_name        => 'AlignSubsetSequence',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::AlignSequence',
      -analysis_capacity => 50,
      -max_retry_count => 0,
      -parameters        => {
        aligner_class  => $aligner_class,
        aligner_dir    => $aligner_dir,
        samtools_dir   => $self->o('samtools_dir'),
        max_intron     => $self->o('max_intron'),
        seq_file_1     => '#subset_seq_file_1#',
        seq_file_2     => '#subset_seq_file_2#',
        store_cmd      => 0,
        threads       => $self->o('threads'),
      },
      -rc_name           => '8GB_multicpu',
      -flow_into         => {
        '2' => {
          'InferStrandness' => {
            'sam_file' => '#sam_file#',
            'sample_seq_file_1' => '#subset_seq_file_1#',
            'sample_seq_file_2' => '#subset_seq_file_2#',
          }
        },
      },
    },

    {
      -logic_name        => 'InferStrandness',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::InferStrandness',
      -parameters        => {
        results_dir    => '#sample_dir#',
        bed_file       => '#genome_bed_file#',
      },
      -rc_name           => 'normal',
      -analysis_capacity => 50,
      -max_retry_count => 0,
      -failed_job_tolerance => 50,
      -flow_into         => {
        '1' => 'AlignSequence',
      },
    },
  );

  my @alignment = (
    {
      -logic_name        => 'AlignSequence',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::AlignSequence',
      -analysis_capacity => 50,
      -max_retry_count => 0,
      -parameters        => {
        aligner_class  => $aligner_class,
        aligner_dir    => $aligner_dir,
        samtools_dir   => $self->o('samtools_dir'),
        max_intron     => $self->o('max_intron'),
        seq_file_1     => '#sample_seq_file_1#',
        seq_file_2     => '#sample_seq_file_2#',
        escape_branch  => -1,
        threads       => $self->o('threads'),
      },
      # Leave at least one retry, so the runnable can pass to the higher runnable
      -max_retry_count => 1,
      -rc_name           => '8GB_multicpu',
      -flow_into         => {
                              '-1' => { 'AlignSequence_HighMem' => INPUT_PLUS() },
                               '2' => { 'SamToBam' => INPUT_PLUS() },
                            },
    },

    {
      -logic_name        => 'AlignSequence_HighMem',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::AlignSequence',
      -analysis_capacity => 30,
      -max_retry_count => 0,
      -can_be_empty      => 1,
      -parameters        => {
        aligner_class  => $aligner_class,
        aligner_dir    => $aligner_dir,
        samtools_dir   => $self->o('samtools_dir'),
        max_intron     => $self->o('max_intron'),
        seq_file_1     => '#sample_seq_file_1#',
        seq_file_2     => '#sample_seq_file_2#',
        threads       => $self->o('threads'),
      },
      -rc_name           => '16GB_multicpu',
      -failed_job_tolerance => 100,
      -flow_into         => {
                               '2' => { 'SamToBam' => INPUT_PLUS() },
                            },
    },

    {
      -logic_name        => 'SamToBam',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::SamToBam',
      -analysis_capacity => 30,
      -max_retry_count => 0,
      -parameters        => {
        samtools_dir   => $self->o('samtools_dir'),
        memory         => $self->o('samtobam_memory'),
        clean_up       => $self->o('clean_up'),
        final_bam_file => '#sample_bam_file#',
        threads       => $self->o('threads'),
      },
      -rc_name           => '16GB_multicpu',
      -flow_into         => {
                              '1->A' => 'Post_alignment',
                              'A->1' => 'WriteBamStats',
                            },
    },
  );

  my @post_alignment = (
    {
      -logic_name => 'Post_alignment',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name    => 'normal',
      -analysis_capacity => 1,
      -flow_into         => WHEN('#dnaseq#', 'DNASeq', ELSE('OrderBam'))
    },

    # DNA-Seq
    {
      -logic_name => 'DNASeq',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name    => 'normal',
      -analysis_capacity => 1,
      -flow_into         => 'DNASeq_BamStats'
    },

    {
      -logic_name        => 'DNASeq_BamStats',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::BamStats',
      -parameters        => {
        bam_file => '#sample_bam_file#',
        results_dir => '#sample_dir#',
      },
      -rc_name           => '8Gb_mem',
      -analysis_capacity => 10,
      -max_retry_count => 0,
      -flow_into         => {
        2 => 'DNASeq_BigWig',
        3 => '?accu_name=bam_stats&accu_address=[]',
      },
    },

    {
      -logic_name        => 'DNASeq_BigWig',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::CreateBigWig',
      -parameters        => {
        bam_file  => '#sample_bam_file#',
      },
      -rc_name           => '16GB',
      -analysis_capacity => 10,
      -max_retry_count => 0,
    },

    # RNA-Seq
    {
      -logic_name        => 'OrderBam',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::OrderBam',
      -parameters        => {
        samtools_dir   => $self->o('samtools_dir'),
        input_bam  => '#sample_bam_file#',
        sorted_bam => '#sorted_bam_file#',
        memory         => 16000,
        threads       => $self->o('threads'),
      },
      -rc_name           => '16GB_multicpu',
      -analysis_capacity => 10,
      -max_retry_count => 0,
      -flow_into         => ['MainBamStats', 'ExtractJunctions']
    },

    {
      -logic_name        => 'MainBamStats',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::BamStats',
      -parameters        => {
        bam_file => '#output_bam_file#',
        results_dir => '#sample_dir#',
      },
      -rc_name           => '8Gb_mem',
      -analysis_capacity => 10,
      -max_retry_count => 0,
      -flow_into         => {
        2 => 'Split_unique',
        3 => '?accu_name=bam_stats&accu_address=[]',
      },
    },

    {
      -logic_name        => 'ExtractJunctions',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::ExtractJunction',
      -parameters => {
        bam_file => '#sorted_bam_file#',
        results_dir => '#sample_dir#',
      },
      -flow_into         => WHEN("-s #genome_gtf_file#", "HtseqFactory"),
      -rc_name           => 'normal',
      -analysis_capacity => 25,
      -max_retry_count => 0,
    },

    {
      -logic_name        => 'HtseqFactory',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::HtseqFactory',
      -parameters        => {
        bam_file => '#sorted_bam_file#',
        is_paired => "#input_is_paired#",
        features => $self->o('features'),
      },
      -rc_name           => 'normal',
      -analysis_capacity => 1,
      -max_retry_count => 0,
      -flow_into         => {
        2 => 'HtseqCount',
      },
    },

    {
      -logic_name        => 'HtseqCount',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::HtseqCountName',
      -parameters        => {
        gtf_file       => '#genome_gtf_file#',
        results_dir    => '#sample_dir#',
      },
      -rc_name           => 'normal',
      -analysis_capacity => 25,
      -max_retry_count => 0,
    },

    {
      -logic_name        => 'Split_unique',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::SplitUnique',
      -parameters => {
        bam_file => '#sample_bam_file#',
      },
      -rc_name           => 'normal',
      -analysis_capacity => 25,
      -max_retry_count => 0,
      -flow_into         => {
        2 => { 'Split_strands' => INPUT_PLUS() },
      },
    },

    {
      -logic_name        => 'Split_strands',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::SplitReadStrand',
      -parameters => {
        bam_file => '#bam_file#',
      },
      -rc_name           => 'normal',
      -analysis_capacity => 25,
      -max_retry_count => 0,
      -flow_into         => {
        2 => 'BamStats',
      },
    },

    {
      -logic_name        => 'BamStats',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::BamStats',
      -parameters        => {
        bam_file => '#bam_file#',
      },
      -rc_name           => '8Gb_mem',
      -analysis_capacity => 10,
      -max_retry_count => 0,
      -flow_into         => {
        1 => { 'CreateBedGraph' => INPUT_PLUS() },
        3 => '?accu_name=bam_stats&accu_address=[]',
      },
    },

    {
      -logic_name        => 'CreateBedGraph',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::CreateBedGraph',
      -parameters        => {
                              bedtools_dir  => $self->o('bedtools_dir'),
                              bam_file => '#sample_bam_file#',
                            },
      -failed_job_tolerance => 100,
      -rc_name           => '8Gb_mem',
      -analysis_capacity => 25,
      -max_retry_count => 0,
      -flow_into         => {
        -1 => 'CreateBedGraph_highmem',
      }
    },

    {
      -logic_name        => 'CreateBedGraph_highmem',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::CreateBedGraph',
      -parameters        => {
                              bedtools_dir  => $self->o('bedtools_dir'),
                              bam_file => '#sample_bam_file#',
                            },
      -failed_job_tolerance => 100,
      -analysis_capacity => 5,
      -max_retry_count => 0,
      -rc_name           => '16Gb_mem',
    },
  );
  
  my @post_process = (

#    {
#      -logic_name => 'Post_process',
#      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
#      -rc_name    => 'normal',
#      -analysis_capacity => 1,
#      -flow_into         => 'WriteBamStats'
#    },
    {
      -logic_name        => 'WriteBamStats',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::WriteBamStats',
      -parameters        => {
        results_dir => '#sample_dir#',
      },
      -rc_name           => 'normal',
      -analysis_capacity => 1,
      -max_retry_count => 0,
      -flow_into         => 'WriteCmdFile',
    },
    {
      -logic_name        => 'WriteCmdFile',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::WriteCmdFile',
      -parameters        => {
                              aligner     => $aligner,
                              results_dir => '#sample_dir#',
                            },
      -rc_name           => 'normal',
      -analysis_capacity => 1,
      -max_retry_count => 0,
      -flow_into         => 'WriteMetadata',
    },
    {
      -logic_name        => 'WriteMetadata',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::WriteMetadata',
      -parameters        => {
                              results_dir => '#sample_dir#',
                            },
      -rc_name           => 'normal',
      -analysis_capacity => 1,
      -max_retry_count => 0,
      -flow_into         => 'Cleanup',
    },
    {
      -logic_name        => 'Cleanup',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::CleanupBam',
      -parameters        => {
                              results_dir => '#sample_dir#',
                              skip_cleanup => $self->o('skip_cleanup'),
                            },
      -rc_name           => 'normal',
      -analysis_capacity => 1,
      -max_retry_count => 0,
      -flow_into         => 'Report',
    },
    {
      -logic_name        => 'Report',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -parameters        => {},
      -rc_name           => 'normal',
      -analysis_capacity => 1,
      -max_retry_count => 0,
    },
  );

  return [
    @backbone,
    @genomes,
    @files,
    @sample_alignment,
    @alignment,
    @post_alignment,
    @post_process
  ];

}

sub resource_classes {
  my ($self) = @_;
  
  my $threads = $self->o("threads");
  
  return {
    %{$self->SUPER::resource_classes},
    '8GB' => {'LSF' => "-q " . $self->o('queue_name') . " -M 8000 -R \"span[hosts=1]\""},
    '16GB' => {'LSF' => "-q " . $self->o('queue_name') . " -M 16000 -R \"span[hosts=1]\""},
    '32GB' => {'LSF' => "-q " . $self->o('queue_name') . " -M 32000 -R \"span[hosts=1]\""},
    '8GB_multicpu' => {'LSF' => "-q " . $self->o('queue_name') . " -n $threads -M 8000 -R \"span[hosts=1]\""},
    '16GB_multicpu' => {'LSF' => "-q " . $self->o('queue_name') . " -n $threads -M 16000 -R \"span[hosts=1]\""},
    '32GB_multicpu' => {'LSF' => "-q " . $self->o('queue_name') . " -n $threads -M 32000 -R \"span[hosts=1]\""},
  }
}

1;
