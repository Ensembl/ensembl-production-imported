## FindPHIBaseCandidates
### [Bio::EnsEMBL::EGPipeline::PipeConfig::FindPHIBaseCandidates_conf](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/FindPHIBaseCandidates_conf.pm)

Adding pathogen-host interactions data from [PHI-base](http://www.phi-base.org) as xrefs to fungal, oomycete and bacterial pathogen genes

### Prerequisites

A registry file with the locations of the core database server(s) and the production database (or `-production_db $PROD_DB_URL` specified).

PHI-base csv file (see below).

### How to run


```
init_pipeline.pl Bio::EnsEMBL::EGPipeline::PipeConfig::FindPHIBaseCandidates_conf \
    $($CMD details script) \
    -hive_force_init 1\
    -registry $REG_FILE \
    -blast_db_dir $BLAST_DB_DIRECTORY \
    -input_file $INPUT_ROTHAMSTED_FILE \
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
| `-input_file` | | path to the prepared  PHI-base current [snapshot](https://github.com/PHI-base/data/blob/master/releases/phi-base_current.csv)
| `-blast_db_dir` | | cache with species-specific (generated from cores) peptide files 
| `-skip_blast_update` | 0 | skip verifying and updating cached BLAST DBs; 1 -- to skip
| `-write` | 0 | remove old and update with new xrefs; 1 -- to actually run (recommended), 0 -- for debugging
| `-phi_release` | | PHI_base release version (i.e. `4.6`)



### Notes

N.B. The pipeline does not support collections at the moment

#### PHI-base data preparation steps

The pipeline works with slightly modified version of [phi-base_current.csv](https://github.com/PHI-base/data/blob/master/releases/phi-base_current.csv) file.

We need a header and unused colums to be removed from the file.

A typical flow for creating such file is:

1. Get the raw [phi-base_current.csv](https://github.com/PHI-base/data/blob/master/releases/phi-base_current.csv) file
```
wget -c https://github.com/PHI-base/data/raw/master/releases/phi-base_current.csv
```

2. Remove headers and unused columns using helper script.
```
# TODO
  > phi_base.filtered.csv
```

3. Normalise ontology terms, using [phi_ontology/phi-base_ontologies.pl](scripts/phi_ontology/phi-base_ontologies.pl) script and ontology maps within [phi_ontology](scripts/phi_ontology) dir.

```
PHIBO_DIR=${ENS_ROOT_DIR}/ensembl-production-imported/scripts/phi_ontology

cp phi_base.filtered.csv phi_base.filtered.csv.orig # creating a back up

perl ${PHIBO_DIR}/phi-base_ontologies.pl \
  -phenotype_dictionary ${PHIBO_DIR}/ontologies_phenotypes.pl \
  -conditions_dictionary ${PHIBO_DIR}/ontologies_conditions.pl \
  -input_file phi_base.filtered.csv \
  -obo_file_destination PHIbase_ontology.${VERSION}.obo
```

Params:
 * `-phenotype_dictionary` -- Perl file with phenotype (`%mapped_phenotypes`) dictionary
 * `-conditions_dictionary` -- Perl file with conditions (`%mapped_conditions`) dictionary
 * `-input_file` -- previous stage result (no header, subset of columns); input file is replaced with the normalised vesrion
 * `-obo_file_destination` --  optional arg to store .obo file

If the script can't solve all the entries in the file,
 it outputs the list of terms 
that need to be added to the dictionaries in [phi_ontology](scripts/phi_ontology) dir.

There are few things to be aware of when running this script:
 * The sript replaces the `input_file` with the normalized version. Please, have a copy if you need one.
 * The script will create (and remove on success) `Ontologies_temp_files` temporary directory,
    in the directory, the `-input_file` is taken from. Thus this recent directory needs to be writable.

4. Now you can use `phi_base.filtered.csv` to initialize `-input_file` option of your pipeline.


### Parts
A few generic from [Common::RunnableDB](docs/Common_RunnableDB.md).

A few from [PHIBase](lib/perl/Bio/EnsEMBL/EGPipeline/PHIBase).

