## RepeatModeler

### Module: [Bio::EnsEMBL::EGPipeline::PipeConfig::RepeatModeler_conf](../lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/RepeatModeler_conf.pm)

Pipeline for building de-novo repeat libraries using [`RepeatModeler`](https://www.repeatmasker.org/RepeatModeler/) utility.

### Prerequisites
A registry file with the locations of the core database server(s).

### How to run

```
# REPEAT_MODELLER_OPTIONS='-max_seq_length 19000000' # or whatever

init_pipeline.pl Bio::EnsEMBL::EGPipeline::PipeConfig::RepeatModeler_conf \
  $($CMD details script) \
  -pipeline_tag  "_${RUN_TAG}" \
  -hive_force_init 1 \
  -registry $REG_FILE \
  -results_dir $RESULTS_DIR \
  -species $SPECIES \
  ${REPEAT_MODELLER_OPTIONS} \
  -do_clustering 1 \
  2> $OUT_DIR/init.stderr \
  1> $OUT_DIR/init.stdout

SYNC_CMD=$(cat $OUT_DIR/init.stdout | grep -- -sync'$' | perl -pe 's/^\s*//; s/"//g')
# should get something like
#   beekeeper.pl -url $url -sync

LOOP_CMD=$(cat $OUT_DIR/init.stdout | grep -- -loop | perl -pe 's/^\s*//; s/\s*#.*$//; s/"//g')
# should get something like
#   beekeeper.pl -url $url -reg_file $REG_FILE -loop

$SYNC_CMD 2> $OUT_DIR/sync.stderr 1> $OUT_DIR/sync.stdout
$LOOP_CMD 2> $OUT_DIR/loop.stderr 1> $OUT_DIR/loop.stdout
```

### Parameters / Options
| option | default value |  meaning | 
| - | - | - |
| `-species` |  | species to process, several `-species ` options are possible
| `-results_dir` | | directory to store results to
| `-min_slice_length` | 5000 | ignore scaffolds that are shorter then this value
| `-max_seq_length` | 10_000_000 | maximum length for a single scaffold, scaffolds will be splitted if exciding this value
| `-max_seqs_per_file` | 10_000 | maximum number of scaffolds to be processed in a chnuck (by a single RepeatModeler) instance
| `-do_clustering` | 0 | deal with redundancy of the final library, using [`CDHit-EST`](http://weizhongli-lab.org/cd-hit/) [clustering](../lib/perl/Bio/EnsEMBL/EGPipeline/DNAFeatures/ClusterRepeatLib.pm)
| `-do_filtering` | 0 | remove repeat models that are found in CDSs ([details](../lib/perl/Bio/EnsEMBL/EGPipeline/DNAFeatures/FilterCustomLib.pm))
| `-repeatmodeler_dir` | | Path to the directory containing the RepeatModeler executables  
| `-blast_engine` | ncbi | which BLAST program to use; options are 'ncbi' or 'wublast' 


### Output
Each chunked portion of the genome will produce a repeat library file.
These are simply concatenated to produce a single file in the `${RESULTS_DIR}` directory,
named `<SPECIES>.rm.lib`, where `<SPECIES>` is the species `production_name`.
There is likely to be some redundancy here, particularly for large genomes that have correspondingly large numbers of chunk files. This isn't ideal, but our main goal is to repeat mask the genome (for which redundancy doesn't matter), so the library is a means to an end, rather than an end in itself. The pipeline uses the multi-core functionality of RepeatModeler (9 cores per hive job), which makes the program faster and enables us to chunk as little as possible; but it remains a necessary evil.
Use `-do_clustering 1` to deal with the library redundancy.

### Parts
A few generic from [Common::RunnableDB](../docs/Common_RunnableDB.md).

A few from [DNAFeatures](../lib/perl/Bio/EnsEMBL/EGPipeline/DNAFeatures/).

