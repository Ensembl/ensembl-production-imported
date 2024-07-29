## SRAAlignment_BRC4 
### Module [Bio::EnsEMBL::EGPipeline::PipeConfig::SRAAlignment_BRC4_conf](../lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/SRAAlignment_BRC4_conf.pm)

Perform RNA(DNA) short read aligments

### Prerequisites

A registry file with the locations of the core database server(s) and the production database (or `-production_db $PROD_DB_URL` specified).

### How to run

```
init_pipeline.pl Bio::EnsEMBL::EGPipeline::PipeConfig::SRAAlignment_BRC4_conf \
    $($CMD details hive) \
    --hive_force_init 1 \
    --reg_conf $REG_FILE \
    --pipeline_dir temp/rnaseq \
    --datasets_file $DATASETS \
    --results_dir $OUTPUT \
    ${OTHER_OPTIONS}
```

### Parameters

| option | default value |  meaning | 
| - | - | - |
| `--pipeline_name` | sra_alignment_brc4_{ensembl_release} | name of the hive pipeline 
| `--results_dir` | | directory where the final alignments are stored
| `--pipeline_dir` | | temp directory
| `--index_dir` | {pipeline_dir}/index | temp directory where the genomes files are extracted
| `--datasets_file` |  | List of datasets following the schema in brc4_rnaseq_schema.json
| `--dna_seq` | 0 | Run the pipeline for DNA-Seq instead of RNA-Seq short reads
| `--global_bw` | 0 | A global BigWig file will be created alongside the analysis
| `--skip_cleanup` |  0, or 1 if using `--dna_seq` | Do not remove intermediate files from the results dir (that includes the bam files)
| `--fallback_ncbi` | 0 | If download from ENA fails, try to download from SRA
| `--force_ncbi` | 0 | Force download data from SRA instead of ENA
| `--sra_dir` | | Where to find fastq-dump (in case it is not in your path)
| `-redo_htseqcount` | 0 | Special pipeline path, where no alignment is done, but a previous alignment is reused to recount features coverage (`--alignments_dir` becomes necessary)
| `--alignments_dir` |  | Where to find previous alignments. Needed for `--redo_htseqcount 1`
| `-infer_metadata` | 1 | Automatically infer the reads metadata: single/paired-end, stranded/unstranded, and in which direction. Otherwise, use the values from the `datasets_file`


### Rare parameters
| option | default value |  meaning | 
| - | - | - |
| `-trim_reads` | 0 | Trim reads before doing any alignment
| `-trimmomatic_bin` | | Trimmomatic path if using `--trim_reads`
| `-trim_adapters_pe` | | Paired-end adapter path if using `--trim_reads`
| `-trim_adapters_se` | | Single-end adapter path if using `--trim_reads`
| `--threads` | 4 | Number of threads to use for various programs
| `--samtobam_memory` | 16,000 | Memory to use for the conversion/sorting SAM -> BAM (in MB)
| `--repeat_masking` | soft | What masking to use when extracting the genomes sequences
| `--max_intron` | 1 | Wether to compute the max intron from the current gene models


### Notes


### Parts
A few generic from [Common::RunnableDB](../docs/Common_RunnableDB.md).

A few from [EnsEMBL::ENA::SRA](../lib/perl/Bio/EnsEMBL/ENA/SRA/).

A few from [BRC4Aligner](../lib/perl/Bio/EnsEMBL/EGPipeline/BRC4Aligner/).

