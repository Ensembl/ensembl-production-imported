## Xref_GPR
### Module [Bio::EnsEMBL::EGPipeline::PipeConfig::Xref_GPR_conf](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/Xref_GPR_conf.pm)

Adding [metabolic reaction annotations from the [Gramene Plant Reactome](https://plantreactome.gramene.org) project to plant genes as xrefs 

### Prerequisites

A registry file with the locations of the core database server(s) and the production database (or `-production_db $PROD_DB_URL` specified).

Reactome flat files (see below).

### How to run

```
init_pipeline.pl Bio::EnsEMBL::EGPipeline::PipeConfig::Xref_GPR_conf \
    $($CMD details script) \
    -hive_force_init 1\
    -registry $REG_FILE \
    -xref_reac_file Ensembl2PlantReactomeReactions.txt \
    -xref_path_file Ensembl2PlantReactome.txt \
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
| `-xref_reac_file` | | flat file with [Gramene Plant Reactome](https://plantreactome.gramene.org) [reactions data](https://plantreactome.gramene.org/download/current/Ensembl2PlantReactomeReactions.txt)
| `-xref_path_file` | | flat file with [Gramene Plant Reactome](https://plantreactome.gramene.org) [pathways data](https://plantreactome.gramene.org/download/current/Ensembl2PlantReactome.txt)
| `-uppercase_gene_id` | 0 | use uppercased versions of stable IDs for mapping; 1 -- to uppercase
| `-production_lookup` | 1 |  Fetch analysis display name, description and web data from the production database; 0 -- to disable


### Notes

The pipeline is designed to be run using
[Ensembl/plant_tools add_gramene_reactome.pl](https://github.com/Ensembl/plant_tools/blob/master/production/core/add_gramene_reactome.pl) script.

Please, see [Ensembl/plant_tools](https://github.com/Ensembl/plant_tools) on how to setup, etc.

Here's the example scenario of usage.

1. Get reactome flat files:
```
wget -c https://plantreactome.gramene.org/download/current/Ensembl2PlantReactomeReactions.txt
wget -c https://plantreactome.gramene.org/download/current/Ensembl2PlantReactome.txt
```

2. Load xrefs using the script itself:
```
ENS_ROOT_DIR=$(pwd) # path to the ensembl repos
ENS_VERSION=104 # or whatever

PIPELINE_DIR=$(pwd)/pipeline_out  # whatever, dir to store intermediate results
mkdir -p ${PIPELINE_DIR}

${ENS_ROOT_DIR}/plant_tools/production/core/add_gramene_reactome.pl
    -v ${ENS_VERSION} \
    -R ${REG_FILE} \
    -P ${PIPELINE_DIR} \
    -reactions Ensembl2PlantReactomeReactions.txt \
    -pathways Ensembl2PlantReactome.txt \
    -s triticum_turgidum \
    -w
```

Sometimes you neef to pass `-uppercase_gene_id 1` option to the underlying pipeline to 
allow usage of uppercase gene stable IDs for mapping (i.e. for _Oryza sativa_ (rice))

3. Make sure you have a proper `analysis_description` entry in you database(s).

### Parts
A few generic from [Common::RunnableDB](docs/Common_RunnableDB.md).

A few from [Xref](lib/perl/Bio/EnsEMBL/EGPipeline/Xref/).

