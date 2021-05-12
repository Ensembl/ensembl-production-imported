## WGA2GenesDirect 
### [Bio::EnsEMBL::EGPipeline::PipeConfig::WGA2GenesDirect_conf](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/WGA2GenesDirect_conf.pm)

Project transripts and create genes based on compara lastz mappings

### Prerequisites

Prepared compara database with LASTZ aligments of the source and target species.

### How to run

```
init_pipeline.pl Bio::EnsEMBL::EGPipeline::PipeConfig::WGA2GenesDirect_conf \
    -hive_force_init 1 \
    -pipeline_tag "_${asm_from}2${asm_to}" \
    $($CMD details script) \
    $($CMD details script_compara_) \
    $($CMD details script_source_) \
    $($CMD details script_target_) \
    $($CMD details script_result_) \
    -compara_dbname $LASTZ_DBNAME \
    -source_dbname $FROM_DBNAME \
    -target_dbname $TO_DBNAME \
    -result_dbname $RES_DBNAME \
    -result_force_rewrite 1 \
    -result_clone_mode 'dna_db' \
    -reg_conf 'none' \
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
| `-pipeline_name` | `wga2genes_${ENS_VERSION}_${pipeline_tag}` | The hive database name will be `${USER}_${pipeline_name}`
| `-pipeline_tag` |  | Tag to append to the  default `-pipeline_name`
| `-source_{host,port,user,pass,dbname}` | | projection source DB connection parameters
| `-target_{host,port,user,pass,dbname}` | | projection target DB connection parameters
| `-compara_{host,port,user,pass,dbname}` | | compara DB with genome alignments (LASTZ) between source and target
| `-result_{host,port,user,pass,dbname}` | | connection details of the DB to store results of the projection

### Notes

N.B. UTRs are usually lost

### Parts

A few from [Bio::EnsEMBL::Analysis](https://github.com/Ensembl/ensembl-analysis/tree/dev/hive_master/modules/Bio/EnsEMBL/Analysis).

