# LICENSE

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

# NAME

SRA Alignment BRC4 conf. A pipeline that aligns reads against genomic
sequences extracted from EnsEMBL core databases.

# SYNOPSIS

    init_pipeline.pl Bio::EnsEMBL::EGPipeline::PipeConfig::SRAAlignment_BRC4_conf \
        --host $HOST --port $PORT --user $USER --pass $PASS \
        --hive_force_init 1 \
        --reg_conf $REG_FILE \
        --pipeline_dir temp/rnaseq \
        --datasets_file $DATASETS \
        --results_dir $OUTPUT \
        ${OTHER_OPTIONS}

Where:

- **--host, --port, --user, --pass**

    Connection details to the MySQL server where the eHive database will be created

- **REG\_FILE**

    An Ensembl registry file, pointing to the core databases to use

- **DATASETS**

    A json representing the datasets to align (following the format from schema/brc4\_rnaseq\_schema.json)

- **OUTPUT**

    The directory where the alignment files will be stored

- **OTHER\_OPTIONS**

    See list of options below

# DESCRIPTION

Perform short read aligments, primarily RNA-Seq (but also supports DNA-Seq).
Cf the ["IN-DEPTH OVERVIEW"](#in-depth-overview) below for more details.

## PARAMETERS

- **--pipeline\_name**

    Name of the hive pipeline.

    Default: sra\_alignment\_brc4\_{ensembl\_release}

- **--results\_dir**

    Directory where the final alignments are stored.

    Mandatory.

- **--pipeline\_dir**

    Temp directory.

    Mandatory.

- **--index\_dir**

    Temp directory where the genomes files are extracted.

    Default: {pipeline\_dir}/index

- **--datasets\_file**

    List of datasets following the schema in brc4\_rnaseq\_schema.json

    Mandatory.

- **--dna\_seq**

    Run the pipeline for DNA-Seq instead of RNA-Seq short reads.

    Default: 0

- **--skip\_cleanup**

    Do not remove intermediate files from the results dir (that includes the bam files).

    0, or 1 if using \`--dna\_seq\`

- **--fallback\_ncbi**

    If download from ENA fails, try to download from SRA.

    Default: 0

- **--force\_ncbi**

    Force download data from SRA instead of ENA.

    Default: 0

- **--sra\_dir**

    Where to find fastq-dump (in case it is not in your path).

- **--redo\_htseqcount**

    Special pipeline path, where no alignment is done, but a previous alignment is reused to recount
    features coverage (\`--alignments\_dir\` becomes necessary).

    Default: 0

- **--alignments\_dir**

    Where to find previous alignments.

    Mandatory if using --redo\_htseqcount

- **--infer\_metadata**

    Automatically infer the reads metadata:

    - single/paired-end
    - stranded/unstranded
    - and in which direction.

    Otherwise, use the values from the \`datasets\_file\`

    Default: 1

## RARE PARAMETERS

- **-trim\_reads**

    Trim reads before doing any alignment.

    Default: 0

- **-trimmomatic\_bin**

    Trimmomatic path if using \`--trim\_reads\`.

- **-trim\_adapters\_pe**

    Paired-end adapter path if using \`--trim\_reads\`.

- **-trim\_adapters\_se**

    Single-end adapter path if using \`--trim\_reads\`

- **--threads**

    Number of threads to use for various programs.

    Default: 4

- **--samtobam\_memory**

    Memory to use for the conversion/sorting SAM -> BAM (in MB).

    Default: 16,000

- **--repeat\_masking**

    What masking to use when extracting the genomes sequences.

    Default: soft

- **--max\_intron**

    Whether to compute the max intron from the current gene models.

    Default: 1

# IN-DEPTH OVERVIEW

- 1. The pipeline first extracts genome data from the cores
- 2. Then it retrieves the runs data from SRA for every sample
- 3. Pre-alignment, each sample data may be trimmed, and their strandness is inferred
- 4. The inferences are checked over the whole of each dataset
- 5. Each sample is then aligned against its reference sequence, and converted into a bam file
- 6. Various post-alignment steps are performed

## 1. GENOME DATA PREPARATION

- Extract the DNA sequence from the core into a fasta file
- Index the fasta file for hisat2 with hisat2-build (without splice sites or exons files)
- Extract the gene models from the core in a gtf format (for htseq-count)
- Extract the gene models from the core in a bed format (for strand inference)

## 2. RNA-SEQ DATA RETRIEVAL

For each run, download the fastq data files from ENA using its SRA accession.

All following processes are performed for each run.

## 3. PRE-ALIGNMENT PROCESSES

### **Trimming**

If the run has a flag 'trim\_reads' set to true, then fastq files are trimmed using trimmomatic
using the following parameters:

    ILLUMINACLIP:$adapters:2:30:10
    LEADING:3
    TRAILING:3
    SLIDINGWINDOW:4:15
    MINLEN:20

Where $adapters is the path to the adapters directory from trimmomatic. There is one directory for
paired-end, and one for single-end reads (the paired/single information must be provided in the
json dataset).

### **Strandness inference**

To ensure that the aligner uses the correct strand information, an additional step to infer the 
strandness of the data is necessary.

Steps:

- 1. Create a subset of reads files with 20,000 reads
- 2. Align those files (without strandness) with hisat2
- 3. Run infer\_experiment.py on the alignment file

The inference compares how the reads are aligned compared to the known gene models:

- If most reads expected to be forward are forward (and vice versa):

        then the data is deemed as stranded in the forward direction.

- If most reads expected to be forward are reversed (and vice versa):

        then the data is deemed as stranded in the reverse direction.

- If the reads are equally in both directions:

        then the data is deemed unstranded.

A cut-off at 85% of aligned reads is applied to discriminate between stranded data:
if the ratio is below this value, then the run is deemed unstranded.

Following hisat2 notation:

- If the data is stranded forward and single-ended, its strandness is stored as "F"
- If the data is stranded reversed and single-ended, its strandness is stored as "R"
- If the data is stranded forward and paired-ended, its strandness is stored as "FR"
- If the data is stranded reversed and paired-ended, its strandness is stored as "RF"

If the strandness or the pair/single-end values differ from those provided in the dataset json file,
then the difference is noted in the log file with a WARNING.
The infer\_experiment output is also stored in this file.

Note that If the values differ, the pipeline continues running using the values inferred.

## 4. ALIGNMENT PARAMETERS CONSENSUS FROM INFERENCES

In order to avoid having datasets with mixed parameters, this step checks all samples inferences
and proposes a consensus.

If a single sample doesn't follow the general consensus, then the runnable Aggregate will fail.

In this case, there are several possibilities:

- Mixed paired-end/single-end

    If this is expected, simply rerun the failing job by adding the following parameter and value:

        ignore_single_paired = 1

- Force a consensus

    If for example there are mixed stranded/unstranded samples in a dataset, you can force all the
    alignments to be unstranded by adding the proposed consensus (output by the failing job), with the
    values replaced, with the following parameter:

        force_aligner_metadata = {consensus}

## 5. ALIGNMENT

Using hisat2 with the following parameters:

    --max-intronlen $max_intron_length # (see below for the value)
    --rna-strandness $strandness # (if stranded, the value is either "F", "R", "FR", or "RF")

The --max-intronlen parameter is computed from the gene set data in the core database,
as the maximum of:

    int(1.5 * $max_intron_length)

where $max\_intron is the longest intron in the gene set

Hisat2 generates an unsorted sam file.

## Temporary Bam files generation

From the file generated by Hisat2, there are two additional conversion steps:

Main bam: sorting by position + conversion to bam file

Name sorted bam: sorting by read name + filter out unaligned reads + conversion to bam file

The first file is temporary and is deleted as the end of the pipeline.

The second file can be conserved and reused for htseq-recount.

The process also generates additional temporary bam files, extracted from the main bam:

- One for unique reads
- One for non-unique reads

If the data is stranded, then each unique/non-unique bam file is also split into:

- Forward stranded reads
- Reverse stranded reads

So if the data is stranded, 4 files are generated. If it is unstranded, 2 files are generated.

## 6. POST-ALIGNMENT PROCESSING

### **Bam stats**

Samtools stats are run on the main bam file, as well as all final split bam files.

The following values are computed:

- coverage (computed with bedtools genomecov)
- mapped (from samtools stats "reads\_mapped" / "raw total sequences")
- number\_reads\_mapped (from samtools stats "reads mapped")
- average\_reads\_length (from samtools stats "average length")
- number\_pairs\_mapped (if paired, from samtools stats "reads properly paired")

Note that for the split bam files, the "mapped" number ratio is over the
main bam "raw total sequences".

### **Create BedGraph**

For each split bam file, the pipeline creates a coverage file in BedGraph format, using:

    bamutils tobedgraph (with stranded parameters --plus or --minus if stranded)
    bedSort

### **Extract junctions**

From the main bam file, the pipeline extracts all the splice junctions into a tabulated text file.

The pipeline extracts the junctions as follow:

- For each read aligned with "N"s in its cigar string
- Get the strand direction from the XS tag
- Get the uniqueness from the NH tag

Create a splice junction for each group of Ns found in the cigar string, with the coordinates,
strand and uniqueness.

### **HTSeq-count**

From the bam file sorted by name, the pipeline runs HTSeq-count:

- Once using unique reads
- Once using all reads

        Note that this file is named "nonunique" because the parameter used
        is "--nonunique all" instead of "--nonunique none")

- If the reads are stranded in the forward direction, use --stranded=yes
- If the reads are stranded in the reverse direction, use --stranded=reverse
- If the reads are not stranded, use --stranded=no
- With parameter --order=name
- With parameter --type=exon --idattr=gene\_id
- With parameter --mode union
- With the gtf file

This is done for each strand, so there will be 4 HTSeq-count files for stranded data,
and 2 for unstranded data.

## Finalisation

Once all files are generated for a run, the pipeline writes two metadata files:

A commands file with the list of commands run (hisat2, etc.)

A metadata file in json that contains in particular the strandness and paired/single end metadata.
This is important for the ulterior htseq-count runs.

At the end of the pipeline, each run contains all its data files in a directory using its run name
defined in the dataset file (note that the name can be any SRA accession like SRX...,
this works as long as there is only one read per accession).

Each run directory contains the following files:

- mappingStats.txt
- metadata.json
- log.txt
- junctions.tab
- commands.json
- 2 or 4 htseq-count files
- 2 or 4 bed files

EBI conserves the final name sorted bam file, but this file is not transferred to UPenn

Each run directory is in a dataset directory.

Each dataset directory is stored in a genome directory using the organism\_abbrev for this genome.

Each genome directory is stored in a component directory (ToxoDB, VectorBase, etc.).

# ALTERNATE ROUTE: DNA-SEQ

It is possible to use this same pipeline to process DNA-Seq datasets.

The main differences are that the post-alignment step only produces a wig file, and keep the
produced bam.

To use this route, simply add the following parameter:

    --dnaseq 1

# ALTERNATE ROUTE: HTSEQ-COUNT RECOUNT

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

- --redo\_htseqcount

    To activate the HTSeq recount route.

- --alignments\_dir

    This is the location of the already aligned datasets (the bam file must be in there, as well as
    the metadata files)  
