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

=head1 NAME

SRA Alignment BRC4 conf. A pipeline that aligns reads against genomic
sequences extracted from EnsEMBL core databases.

=head1 SYNOPSIS

  init_pipeline.pl Bio::EnsEMBL::EGPipeline::PipeConfig::SRAAlignment_BRC4_conf \
      --host $HOST --port $PORT --user $USER --pass $PASS \
      --hive_force_init 1 \
      --reg_conf $REG_FILE \
      --pipeline_dir temp/rnaseq \
      --datasets_file $DATASETS \
      --results_dir $OUTPUT \
      ${OTHER_OPTIONS}

Where:

=over

=item B<--host, --port, --user, --pass>

Connection details to the MySQL server where the eHive database will be created

=item B<REG_FILE>

An Ensembl registry file, pointing to the core databases to use

=item B<DATASETS>

A json representing the datasets to align (following the format from schema/brc4_rnaseq_schema.json)

=item B<OUTPUT>

The directory where the alignment files will be stored

=item B<OTHER_OPTIONS>

See list of options below

=back

=head1 DESCRIPTION

Perform short read aligments, primarily RNA-Seq (but also supports DNA-Seq).
Cf the L<"IN-DEPTH OVERVIEW"> below for more details.

=head2 PARAMETERS

=over

=item B<--pipeline_name>

Name of the hive pipeline.

Default: sra_alignment_brc4_{ensembl_release}

=item B<--results_dir>

Directory where the final alignments are stored.

Mandatory.

=item B<--pipeline_dir>

Temp directory.

Mandatory.

=item B<--index_dir>

Temp directory where the genomes files are extracted.

Default: {pipeline_dir}/index

=item B<--datasets_file>

List of datasets following the schema in brc4_rnaseq_schema.json

Mandatory.

=item B<--dna_seq>

Run the pipeline for DNA-Seq instead of RNA-Seq short reads.

Default: 0

=item B<--skip_cleanup>

Do not remove intermediate files from the results dir (that includes the bam files).

0, or 1 if using `--dna_seq`

=item B<--fallback_ncbi>

If download from ENA fails, try to download from SRA.

Default: 0

=item B<--force_ncbi>

Force download data from SRA instead of ENA.

Default: 0

=item B<--sra_dir>

Where to find fastq-dump (in case it is not in your path).

=item B<--redo_htseqcount>

Special pipeline path, where no alignment is done, but a previous alignment is reused to recount
features coverage (`--alignments_dir` becomes necessary).

Default: 0

=item B<--alignments_dir>

Where to find previous alignments.

Mandatory if using --redo_htseqcount

=item B<--infer_metadata>

Automatically infer the reads metadata:

=over

=item * single/paired-end

=item * stranded/unstranded

=item * and in which direction.

=back

Otherwise, use the values from the `datasets_file`

Default: 1

=back

=head2 RARE PARAMETERS

=over

=item B<-trim_reads>

Trim reads before doing any alignment.

Default: 0

=item B<-trimmomatic_bin>

Trimmomatic path if using `--trim_reads`.

=item B<-trim_adapters_pe>

Paired-end adapter path if using `--trim_reads`.

=item B<-trim_adapters_se>

Single-end adapter path if using `--trim_reads`

=item B<--threads>

Number of threads to use for various programs.

Default: 4

=item B<--samtobam_memory>

Memory to use for the conversion/sorting SAM -> BAM (in MB).

Default: 16,000

=item B<--repeat_masking>

What masking to use when extracting the genomes sequences.

Default: soft

=item B<--max_intron>

Whether to compute the max intron from the current gene models.

Default: 1

=back

=head1 IN-DEPTH OVERVIEW

=over

=item 1) The pipeline first extracts genome data from the cores

=item 2) Then it retrieves the runs data from SRA for every sample

=item 3) Pre-alignment, each sample data may be trimmed, and their strandness is inferred

=item 4) The inferences are checked over the whole of each dataset

=item 5) Each sample is then aligned against its reference sequence, and converted into a bam file

=item 6) Various post-alignment steps are performed

=back

=head2 1) GENOME DATA PREPARATION

=over

=item * Extract the DNA sequence from the core into a fasta file

=item * Index the fasta file for hisat2 with hisat2-build (without splice sites or exons files)

=item * Extract the gene models from the core in a gtf format (for htseq-count)

=item * Extract the gene models from the core in a bed format (for strand inference)

=back

=head2 2) RNA-SEQ DATA RETRIEVAL

For each run, download the fastq data files from ENA using its SRA accession.

All following processes are performed for each run.

=head2 3) PRE-ALIGNMENT PROCESSES

=head3 B<Trimming>

If the run has a flag 'trim_reads' set to true, then fastq files are trimmed using trimmomatic
using the following parameters:

    ILLUMINACLIP:$adapters:2:30:10
    LEADING:3
    TRAILING:3
    SLIDINGWINDOW:4:15
    MINLEN:20

Where $adapters is the path to the adapters directory from trimmomatic. There is one directory for
paired-end, and one for single-end reads (the paired/single information must be provided in the
json dataset).

=head3 B<Strandness inference>

To ensure that the aligner uses the correct strand information, an additional step to infer the 
strandness of the data is necessary.

Steps:

=over

=item 1) Create a subset of reads files with 20,000 reads

=item 2) Align those files (without strandness) with hisat2

=item 3) Run infer_experiment.py on the alignment file

=back

The inference compares how the reads are aligned compared to the known gene models:

=over

=item * If most reads expected to be forward are forward (and vice versa):

then the data is deemed as stranded in the forward direction.

=item * If most reads expected to be forward are reversed (and vice versa):

then the data is deemed as stranded in the reverse direction.

=item * If the reads are equally in both directions:

then the data is deemed unstranded.

=back

A cut-off at 85% of aligned reads is applied to discriminate between stranded data:
if the ratio is below this value, then the run is deemed unstranded.


Following hisat2 notation:

=over

=item * If the data is stranded forward and single-ended, its strandness is stored as "F"

=item * If the data is stranded reversed and single-ended, its strandness is stored as "R"

=item * If the data is stranded forward and paired-ended, its strandness is stored as "FR"

=item * If the data is stranded reversed and paired-ended, its strandness is stored as "RF"

=back

If the strandness or the pair/single-end values differ from those provided in the dataset json file,
then the difference is noted in the log file with a WARNING.
The infer_experiment output is also stored in this file.

Note that If the values differ, the pipeline continues running using the values inferred.

=head2 4) ALIGNMENT PARAMETERS CONSENSUS FROM INFERENCES

In order to avoid having datasets with mixed parameters, this step checks all samples inferences
and proposes a consensus.

If a single sample doesn't follow the general consensus, then the runnable Aggregate will fail.

In this case, there are several possibilities:

=over

=item Mixed paired-end/single-end

If this is expected, simply rerun the failing job by adding the following parameter and value:

  ignore_single_paired = 1

=item Force a consensus

If for example there are mixed stranded/unstranded samples in a dataset, you can force all the
alignments to be unstranded by adding the proposed consensus (output by the failing job), with the
values replaced, with the following parameter:

  force_aligner_metadata = {consensus}

=back

=head2 5) ALIGNMENT

Using hisat2 with the following parameters:

  --max-intronlen $max_intron_length # (see below for the value)
  --rna-strandness $strandness # (if stranded, the value is either "F", "R", "FR", or "RF")

The --max-intronlen parameter is computed from the gene set data in the core database,
as the maximum of:

  int(1.5 * $max_intron_length)

where $max_intron is the longest intron in the gene set

Hisat2 generates an unsorted sam file.

=head2 Temporary Bam files generation

From the file generated by Hisat2, there are two additional conversion steps:

Main bam: sorting by position + conversion to bam file

Name sorted bam: sorting by read name + filter out unaligned reads + conversion to bam file

The first file is temporary and is deleted as the end of the pipeline.

The second file can be conserved and reused for htseq-recount.

The process also generates additional temporary bam files, extracted from the main bam:

=over

=item One for unique reads

=item One for non-unique reads

=back

If the data is stranded, then each unique/non-unique bam file is also split into:

=over

=item Forward stranded reads

=item Reverse stranded reads

=back

So if the data is stranded, 4 files are generated. If it is unstranded, 2 files are generated.

=head2 6) POST-ALIGNMENT PROCESSING

=head3 B<Bam stats>

Samtools stats are run on the main bam file, as well as all final split bam files.

The following values are computed:

=over

=item * coverage (computed with bedtools genomecov)

=item * mapped (from samtools stats "reads_mapped" / "raw total sequences")

=item * number_reads_mapped (from samtools stats "reads mapped")

=item * average_reads_length (from samtools stats "average length")

=item * number_pairs_mapped (if paired, from samtools stats "reads properly paired")

=back

Note that for the split bam files, the "mapped" number ratio is over the
main bam "raw total sequences".

=head3 B<Create BedGraph>

For each split bam file, the pipeline creates a coverage file in BedGraph format, using:

  bamutils tobedgraph # (with stranded parameters --plus or --minus if stranded)
  bedSort

=head3 B<Extract junctions>

From the main bam file, the pipeline extracts all the splice junctions into a tabulated text file.

The pipeline extracts the junctions as follow:

=over

=item * For each read aligned with "N"s in its cigar string

=item * Get the strand direction from the XS tag

=item * Get the uniqueness from the NH tag

=back

Create a splice junction for each group of Ns found in the cigar string, with the coordinates,
strand and uniqueness.

=head3 B<HTSeq-count>

From the bam file sorted by name, the pipeline runs HTSeq-count:

=over

=item * Once using unique reads

=item * Once using all reads

Note that this file is named "nonunique" because the parameter used is "--nonunique all"
instead of "--nonunique none")

=item * If the reads are stranded in the forward direction, use --stranded=yes

=item * If the reads are stranded in the reverse direction, use --stranded=reverse

=item * If the reads are not stranded, use --stranded=no

=item * With parameter --order=name

=item * With parameter --type=exon --idattr=gene_id

=item * With parameter --mode union

=item * With the gtf file

=back

This is done for each strand, so there will be 4 HTSeq-count files for stranded data,
and 2 for unstranded data.

=head2 Finalisation

Once all files are generated for a run, the pipeline writes two metadata files:

A commands file with the list of commands run (hisat2, etc.)

A metadata file in json that contains in particular the strandness and paired/single end metadata.
This is important for the ulterior htseq-count runs.


At the end of the pipeline, each run contains all its data files in a directory using its run name
defined in the dataset file (note that the name can be any SRA accession like SRX...,
this works as long as there is only one read per accession).

Each run directory contains the following files:

=over

=item * mappingStats.txt

=item * metadata.json

=item * log.txt

=item * junctions.tab

=item * commands.json

=item * 2 or 4 htseq-count files

=item * 2 or 4 bed files

=back

EBI conserves the final name sorted bam file, but this file is not transferred to UPenn

Each run directory is in a dataset directory.

Each dataset directory is stored in a genome directory using the organism_abbrev for this genome.

Each genome directory is stored in a component directory (ToxoDB, VectorBase, etc.).

=head1 ALTERNATE ROUTE: DNA-SEQ

It is possible to use this same pipeline to process DNA-Seq datasets.

The main differences are that the post-alignment step only produces a wig file, and keep the
produced bam.

To use this route, simply add the following parameter:

  --dnaseq 1

=head1 ALTERNATE ROUTE: HTSEQ-COUNT RECOUNT

If you have already aligned a dataset, and the corresponding genome has had a gene set update, then
you can rerun only the recount part of the pipeline. To do so, use the following parameters:

  init_pipeline.pl Bio::EnsEMBL::EGPipeline::PipeConfig::SRAAlignment_BRC4_conf \
  --host $HOST --port $PORT --user $USER --pass $PASS \
  --pipeline_name rnaseq_recount \
  --pipeline_dir temp/rnaseq_recount \
  --reg_conf $REGISTRY \
  --results_dir $OUTPUT \
  --datasets_file $DATASETS \
  --redo_htseqcount 1 \
  --alignments_dir $ALIGNMENTS_DIR

Note the addition of 2 parameters:

=over

=item --redo_htseqcount

To activate the HTSeq recount route.

=item --alignments_dir

This is the location of the already aligned datasets (the bam file must be in there, as well as
the metadata files)  

=back


=cut

package Bio::EnsEMBL::EGPipeline::PipeConfig::SRAAlignment_BRC4_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf');

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use Bio::EnsEMBL::Hive::Version 2.5;

use File::Spec::Functions qw(catdir catfile);

use File::Basename;
use Class::Inspector;
my $package_path = Class::Inspector->loaded_filename(__PACKAGE__);
my $package_dir = dirname($package_path);

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
    datasets_json_schema => catfile($package_dir, '../BRC4Aligner/brc4_rnaseq_schema.json'),
    
    # BEHAVIOUR
    # Use DNA-Seq specific parameters
    # In that case, do not remove bam files
    dnaseq => 0,
    skip_cleanup => $self->o('dnaseq'),
    
    # Download via SRA from NCBI if ENA fails (fallback)
    fallback_ncbi => 0,
    sra_dir => undef,
    # Force using NCBI SRA instead of ENA
    force_ncbi => 0,
    
    # If there is already an alignment for a sample, redo its htseq-count only
    # (do nothing otherwise)
    redo_htseqcount => 0,
    features => ['exon'],
    alignments_dir => undef,
    
    # Use input metadata, instead of inferring them (pair/strand)
    infer_metadata => 1,

    # For heavy analyses, use multiple cpus
    threads    => 4,
    samtobam_memory => 16000,

    # Do not proceed if the reads are too long
    max_read_length => 1000,

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

sub hive_meta_table {
    my ($self) = @_;
    return {
        %{$self->SUPER::hive_meta_table},       # here we inherit anything from the base class
        'hive_use_param_stack'  => 1,           # switch on the new param_stack mechanism
    };
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
    'force_ncbi'      => $self->o('force_ncbi'),
  };
}

sub pipeline_analyses {
  my ($self) = @_;

  my $aligner = 'hisat2';
  my $aligner_class = 'Bio::EnsEMBL::EGPipeline::Common::Aligner::HISAT2Aligner';
  my $aligner_dir = $self->o('hisat2_dir');

  my @backbone = (
    {
      -logic_name        => 'Init',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'mkdir -p ' . $self->o('pipeline_dir') . "; " .
               'mkdir -p ' . $self->o('results_dir'),
      },
      -input_ids  => [{}],
      -rc_name           => 'normal',
      -meadow_type       => 'LSF',
      -analysis_capacity => 1,
      -max_retry_count => 0,
      -flow_into  => { 1 => 'Check_schema' },
    },


    {
      -logic_name => 'Check_schema',
      -module     => 'ensembl.pipeline.json.schema_validator',
      -language     => 'python3',
      -parameters => {
        json_file => $self->o('datasets_file'),
        json_schema => $self->o('datasets_json_schema'),
      },
      -flow_into  => {
        '1->A' => 'Dataset_species_factory',
        'A->1' => 'Email_report',
      },
      -max_retry_count => 0,
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
      -analysis_capacity => 1,
    },

    {
      -logic_name => 'Dataset_species_factory',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::DatasetSpeciesFactory',
      -parameters => {
        datasets_file => $self->o('datasets_file'),
      },
      -flow_into  => {
        '2->A' => 'Prepare_genome',
        'A->2' => 'Species_report',
        '3' => 'Organisms_not_found',
      },
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
      -analysis_capacity => 1,
      -max_retry_count => 0,
    },

    {
      -logic_name        => 'Organisms_not_found',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name           => 'normal',
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
        1 => 'Genome_data_factory',
      },
    },

    {
      -logic_name        => 'Genome_data_factory',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name           => 'normal',
      -flow_into         => {
        '1->A' => 'Genome_factory',
        'A->1' => 'Datasets_factory',
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
        WHEN('not -f #genome_file# . ".indexed"', 'Dump_and_index'),
      ],
    },

    {
      -logic_name        => 'Dump_and_index',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name           => 'normal',
      -flow_into         => {
        '1->A' => 'Dump_genome',
        'A->1' => 'Finished_dump',
      },
    },

    {
      -logic_name        => 'Finished_dump',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'touch #genome_file#.indexed',
      },
      -rc_name           => 'normal',
      -max_retry_count => 0,
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
    {
      -logic_name        => 'Datasets_factory',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::DatasetFactory',
      -parameters        => {
        results_dir => '#species_results_dir#',
        datasets_file => $self->o('datasets_file'),
        tax_id_restrict => 0,
        check_library => 0,
      },
      -failed_job_tolerance => 0,
      -analysis_capacity => 1,
      -max_retry_count => 0,
      -rc_name           => 'normal',
      -flow_into         => {
        '2->A' => WHEN("not #redo_htseqcount#", 'Sub_Samples_factory'),
        'A->2' => WHEN("not #redo_htseqcount#", 'Aggregate_metadata'),
        '2' => WHEN("#redo_htseqcount#", 'Redo_Samples_factory'),
        '3' => 'Datasets_not_redone',
      },
    },

    {
      -logic_name        => 'Datasets_not_redone',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name           => 'normal',
      -analysis_capacity => 1,
      -max_retry_count => 0,
    },

    {
      -logic_name        => 'Aggregate_metadata',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::AggregateMetadata',
      -failed_job_tolerance => 100,
      -analysis_capacity => 1,
      -max_retry_count => 0,
      -rc_name           => 'normal',
      -flow_into         => {
        '2' => 'Samples_factory'
      },
    },

    {
      -logic_name        => 'Samples_factory',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::SampleFactory',
      -parameters        => {
        results_dir => '#species_results_dir#',
        tax_id_restrict => 0,
        check_library => 0,
      },
      -failed_job_tolerance => 0,
      -analysis_capacity => 1,
      -max_retry_count => 0,
      -rc_name           => 'normal',
      -flow_into         => {
        '4' => 'Prepare_fastq'
      },
    },

    {
      -logic_name        => 'Prepare_fastq',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::PrepareFastq',
      -parameters        => {
        work_dir => '#species_work_dir#',
      },
      -failed_job_tolerance => 10,
      -analysis_capacity => 1,
      -max_retry_count => 0,
      -rc_name           => 'normal',
      -flow_into         => {
        2 => "AlignSequence",
      }
    },

    {
      -logic_name        => 'Sub_Samples_factory',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::SampleFactory',
      -parameters        => {
        results_dir => '#species_results_dir#',
        tax_id_restrict => 0,
        check_library => 0,
      },
      -failed_job_tolerance => 0,
      -analysis_capacity => 1,
      -max_retry_count => 0,
      -rc_name           => 'normal',
      -flow_into         => {
        '2->A' => 'Download',
        'A->4' => 'Sub_Merge_fastq',
      },
    },

    {
      -logic_name        => 'Download',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name           => 'normal',
      -analysis_capacity => 1,
      -flow_into         => {
        '1' => WHEN('#force_ncbi#', 'SRASeqFileFromNCBI', ELSE 'SRASeqFileFromENA'),
      },
    },

    # HTseq recount path
    {
      -logic_name        => 'Redo_Samples_factory',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::SampleFactory',
      -parameters        => {
        results_dir => '#species_results_dir#',
        tax_id_restrict => 0,
        check_library => 0,
      },
      -failed_job_tolerance => 0,
      -analysis_capacity => 1,
      -max_retry_count => 0,
      -rc_name           => 'normal',
      -flow_into         => {
        '3' => 'GenesCheckForRedo',
      },
    },

    {
      -logic_name        => 'GenesCheckForRedo',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::GenesCheck',
      -analysis_capacity => 1,
      -max_retry_count => 0,
      -rc_name           => 'normal',
      -flow_into         => {
        2 => "SyncAlignmentFiles",
      },
    },

    {
      -logic_name        => 'SyncAlignmentFiles',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::SyncAlignmentFiles',
      -analysis_capacity => 2,
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
        2 => 'HtseqFactory_redo',
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
        2 => 'HtseqCount_redo',
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
      -logic_name        => 'Sub_Merge_fastq',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::MergeFastq',
      -parameters        => {
        work_dir => '#species_work_dir#',
      },
      -failed_job_tolerance => 10,
      -analysis_capacity => 1,
      -max_retry_count => 0,
      -rc_name           => 'normal',
      -flow_into         => {
        2 => "PreAlignment",
      }
    },

    {
      -logic_name        => 'SRASeqFileFromENA',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::SRASeqFile',
      -failed_job_tolerance => 10,
      -analysis_capacity => 4,
      -max_retry_count => 3,
      -failed_job_tolerance => 10,
      -parameters        => {
        work_dir => '#species_work_dir#',
        fallback_ncbi => $self->o('fallback_ncbi'),
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
        max_read_length => $self->o('max_read_length'),
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
        '1' => WHEN('#trim_reads#', 'Trim', ELSE 'GetMetadata')
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
        '2' => 'GetMetadata',
      },
    },
  );

    ####################################################################
    # Get metadata for each sample, to use for the aligner
  my @get_metadata = (

    {
      -logic_name        => 'GetMetadata',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name           => 'normal',
      -analysis_capacity => 1,
      -max_retry_count => 0,
      -flow_into         => {
        '1->A' => 'CheckMetadata',
        'A->1' => 'SendMetadata',
      },
    },
    {
      -logic_name        => 'SendMetadata',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name           => 'normal',
      -analysis_capacity => 1,
      -max_retry_count => 0,
      -flow_into         => {
        1 => [
          '?accu_name=aligner_metadata_hash&accu_input_variable=aligner_metadata&accu_address={sample_name}',
        ],
      },
    },


    {
      -logic_name        => 'CheckMetadata',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::GenesCheck',
      -analysis_capacity => 1,
      -max_retry_count => 0,
      -rc_name           => 'normal',
      -flow_into         => {
        2 => WHEN("#infer_metadata#", "SubsetSequence", ELSE("UseInputMetadata")),
        3 => "UseInputMetadata",
      },
    },

    {
      -logic_name        => 'UseInputMetadata',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name           => 'normal',
      -analysis_capacity => 1,
      -max_retry_count => 0,
      -flow_into         => {
        1 => '?accu_name=aligner_metadata&accu_input_variable=input_metadata',
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
        '1' => 'AlignSubsetSequence',
      },
    },

    {
      -logic_name        => 'AlignSubsetSequence',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::AlignSequence',
      -failed_job_tolerance => 10,
      -analysis_capacity => 50,
      -max_retry_count => 0,
      -parameters        => {
        aligner_class  => $aligner_class,
        aligner_dir    => $aligner_dir,
        samtools_dir   => $self->o('samtools_dir'),
        max_intron     => $self->o('max_intron'),
        seq_file_1     => '#subset_seq_file_1#',
        seq_file_2     => '#subset_seq_file_2#',
        aligner_metadata => '#input_metadata#',
        store_cmd      => 0,
        threads       => 1,
      },
      -rc_name           => '8GB',
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
        2 => '?accu_name=aligner_metadata&accu_input_variable=aligner_metadata',
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
        aligner_metadata => '#aggregated_aligner_metadata#',
        escape_branch  => -1,
        threads       => $self->o('threads'),
      },
      # Leave at least one retry, so the runnable can pass to the higher runnable
      -max_retry_count => 1,
      -rc_name           => '8GB_multicpu',
      -flow_into         => {
                              '-1' => 'AlignSequence_HighMem',
                               '2' => 'SamToBam',
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
        aligner_metadata => '#aggregated_aligner_metadata#',
        threads       => $self->o('threads'),
      },
      -rc_name           => '16GB_multicpu',
      -failed_job_tolerance => 100,
      -flow_into         => {
                               '2' => 'SamToBam',
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
        skip_cleanup   => $self->o('skip_cleanup'),
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
        aligner_metadata => '#aggregated_aligner_metadata#',
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
        aligner_metadata => '#aggregated_aligner_metadata#',
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
        aligner_metadata => '#aggregated_aligner_metadata#',
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
        aligner_metadata => '#aggregated_aligner_metadata#',
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
        aligner_metadata => '#aggregated_aligner_metadata#',
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
        aligner_metadata => '#aggregated_aligner_metadata#',
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
        2 => 'Split_strands',
      },
    },

    {
      -logic_name        => 'Split_strands',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::SplitReadStrand',
      -parameters => {
        aligner_metadata => '#aggregated_aligner_metadata#',
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
        aligner_metadata => '#aggregated_aligner_metadata#',
        bam_file => '#bam_file#',
      },
      -rc_name           => '8Gb_mem',
      -analysis_capacity => 10,
      -max_retry_count => 0,
      -flow_into         => {
        1 => 'CreateBedGraph',
        3 => '?accu_name=bam_stats&accu_address=[]',
      },
    },

    {
      -logic_name        => 'CreateBedGraph',
      -module            => 'Bio::EnsEMBL::EGPipeline::BRC4Aligner::CreateBedGraph',
      -parameters        => {
        aligner_metadata => '#aggregated_aligner_metadata#',
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
        aligner_metadata => '#aggregated_aligner_metadata#',
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
    @get_metadata,
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
