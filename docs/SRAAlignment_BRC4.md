## SRAAlignment_BRC4 
### Module [Bio::EnsEMBL::EGPipeline::PipeConfig::SRAAlignment_BRC4_conf](../lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/SRAAlignment_BRC4_conf.pm)

Perform RNA(DNA) short read aligments

### Prerequisites

A registry file with the locations of the core database server(s) and the production database (or `-production_db $PROD_DB_URL` specified).

### How to run

```
init_pipeline.pl Bio::EnsEMBL::EGPipeline::PipeConfig::SRAAlignment_BRC4_conf \
    $($CMD details script) \
    -hive_force_init 1\
    -registry $REG_FILE \
    ${OTHER_OPTIONS} \
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
| `-pipeline_dir` | | directory to store results to
| `-meta_filters` | |
| `-dnaseq` | 0 | 
| `-use_ncbi` | 0 | 
| `-clean_up` | 0 | 
| `-redo_htseqcount` | 0 |
| `-infer_metadata` | 1 |
| `-trimmomatic_bin` |
| `-trim_adapters_pe` |
| `-trim_adapters_se` |
| `-repeat_masking` | `soft` | 

### Notes


### Parts
A few generic from [Common::RunnableDB](../docs/Common_RunnableDB.md).

A few from [EnsEMBL::ENA::SRA](../lib/perl/Bio/EnsEMBL/ENA/SRA/).

A few from [BRC4Aligner](../lib/perl/Bio/EnsEMBL/EGPipeline/BRC4Aligner/).

