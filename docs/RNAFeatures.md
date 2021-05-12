## RNAFeatures
### Module [Bio::EnsEMBL::EGPipeline::PipeConfig::RNAFeatures_conf](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/RNAFeatures_conf.pm)

The RNA features pipeline aligns RNA covariance models against a genome.
By default it runs cmscan with the latest Rfam covariance models, tRNAscan for tRNA, and loads miRNA genes from miRBase.
If you have other sets of covariance models, they can also be passed to cmscan. The results are stored as alignments

### Prerequisites

A registry file with the locations of the core database server(s) and the production database (or `-production_db $PROD_DB_URL` specified).

Prepared Rfam data (see below).

### How to run

```
init_pipeline.pl Bio::EnsEMBL::EGPipeline::PipeConfig::RNAFeatures_conf \
    $($CMD details script) \
    -hive_force_init 1\
    -registry $REG_FILE \
    -production_db "$($PROD_SERVER details url)""$PROD_DBNAME" \
    -pipeline_dir $OUT_DIR/rna_features \
    -species $SPECIES \
    -eg_pipelines_dir $ENS_DIR/ensembl-production-imported \
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
| `-run_trnascan` | 1 | run `tRNAscan`; 0 -- to disable
| `-run_cmscan` | 1 | run `cmscan`
| `-rfam_rrna` | 1 | include rRNA set of Rfam covariance models; 0 -- to exclude 
| `-rfam_trna` | 0 | include tRNA set of Rfam covariance models; 1 -- to include 
| `-cmscan_cm_file ${SPECIES}\|_all_=/path/to/lib/file` | | Custom set of models to use instead of default Rfam covariance models (`-rfam_cm_file`); can be used several time; _all_ to use the same for all species
| `-cmscan_logic_name ${SPECIES}\|_all_=` | all=cmscan_custom |  Set custom logic name; default can be overidden on a per-species basis (several options), or by using `all=` -- for all species
| `-cmscan_db_name ${SPECIES}\|_all_=${ALT_DB_NAME}` | `all=Rfam` | use `external_db` (pre-existinf in prod db) to associate alignments with
| `-cmscan_cpu` | 3 | Number of cores to be use by `cmscan` (sets LSF's `-n` option)
| `-cmscan_heuristics` | default | sets trade-off between sensitivity and run-time; can be one of `slowest`, `slower`, `slow`, `default`, `faster`, `fastest`, `rfam_official` (recommended for bacteria runs) (see [Bio::EnsEMBL::Analysis::Runnable::CMScan](lib/perl/Bio/EnsEMBL/Analysis/Runnable/CMScan.pm) and [cmscan man](https://manpages.ubuntu.com/manpages/xenial/man1/cmscan.1.html) for mappings
| `-cmscan_threshold` | 0.001 | E-value threshold for reporting alignments
| `-taxonomic_filtering` | 1 | Enable taxonomic filtering of Rfam models (see below for details); 0 -- to disable
| `-taxonomic_threshold` |  0.02 | (fraction 0 to 1; effective when `-taxonomic_filtering 1`) Only allow Rfam model if ratio of number of species within the devision taxonomic level (or `-taxonomic_levels` parameter) to the total number of species associated with this model is higher than threshold
| `-rfam_taxonomy_file` | `${rfam_dir}/taxonomic_levels.txt` | Path to the file with counts of Rfam sequences at division-specific taxonomic levels
| `-taxonomic_levels` | | List of taxomomic levels to be used for taxonomic filtering 
| `-taxonomic_lca` | 0 | Enable LCA in addition to the division level filtering; 1 -- to enable 
| `-load_mirbase` | 1 | look for miRBase files, if there are none then nothing will happen; 0 -- to disable 
| `-mirbase_file ${SPECIES}=/path/to/file.gff3` | | specify which [miRBase](https://www.mirbase.org/) file to use for which species; should be loaded in advance (see mirBase [species range](ftp://mirbase.org/pub/mirbase/CURRENT/genomes/)); can be outdated, check assembly versions 
| `-delete_existing` | 1 | Delete pre-existing analysis data; 0 to disable
| `-max_hive_capacity` | 50 | The default hive capacity (i.e. the maximum number of concurrent workers) (no enforced upper limit, but above 150 you might run into problems with locking or database connections).
| `-cmscan_exe` | `cmscan` |  Path to the cmscan executable file.
| `-cmscan_parameters` | | Parameters that are passed directly to the cmscan executable
| `-rfam_version` | see [RNAFeatures_conf.pm](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/RNAFeatures_conf.pm) | set `rfam_version` to be included into default `-rfam_dir` value (see [RNAFeatures_conf.pm](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/RNAFeatures_conf.pm))
| `-rfam_dir` | [RFAM_VERSIONS_DIR](lib/perl/Bio/EnsEMBL/EGPipeline/PrivateConfDetails/Impl.pm.example)/`${rfam_version}` | path to directory containing Rfam files
| `-rfam_cm_file` | `${rfam_dir}/Rfam.cm` | Path to the Rfam CM file
| `-rfam_logic_name` | `cmscan_rfam_${rfam_version}` | Logic name for the Rfam analysis
| `-rfam_db_name` | `RFAM` | external_db name of the Rfam record
| `-trnascan_exe` | `tRNAscan-SE` | Path to the trnascan executable file
| `-trnascan_logic_name` | `trnascan_align` | logic name for the the tRNAscan analysis
| `-trnascan_db_name` | `TRNASCAN_SE` | external_db name of the tRNAscan 
| `-trnascan_pseudo` | 0 | include prediction of tRNA pseudogenes; 1 -- to include
| `-trnascan_parameters` | | parameters passed to the trnascan executable
| `-mirbase_logic_name` | `mirbase` | logic name for the the miRBase analysis
| `-mirbase_db_name` | `miRBase` | external_db name of the miRBase
| `-mirbase_version` | see [RNAFeatures_conf.pm](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/RNAFeatures_conf.pm) | version of miRBase files
| `-pipeline_name` | `rna_features_${ENS_VERSION}` | The hive database name will be `${USER}_${pipeline_name}`
| `-production_lookup` | 1 |  Fetch analysis display name, description and web data from the production database; 0 -- to disable


### Notes

#### Taxonomic filtering

Most of the covariance models (CMs) in Rfam have been constructed with reference to a particular set of species.
Thus allowing primates only CMs not to be annotated in protists.
But in practice some RNA genes (being derived from transposons) get many hundreds of hits in inappropriate species.

Rfam doesn't provide information on the taxonomic restrictedness of CMs,
but the taxonomic spread of the sequences associated with each model can be calculated.
This data has been precalculated at the division-level (the default) level of filtering.

To decide if an Rfam model can be applied to a particular species,
the pipeline checks whether >2% (see `-taxonomic_threshold` parameter) of the species associated with the model
are within the taxonomic level equivalent to the species' division.

Taxonomic filtering is switched on by default, but can be switched off with the `-taxonomic_filtering` parameter.
The 2% threshold can be changed via the `-taxonomic_threshold` parameter, which takes values from 0 to 1.

#### Strict taxonomic filtering

The pipeline allows you to filter on LCA in addition to the division level filtering,
using the `-taxonomic_lca parameter` (default 0).
This level of filtering is in addition to the threshold mentioned above, and the blacklists and whitelists mentioned below.

Using the LCA option appends `_lca` to the analysis `logic_name`.
Thus the pipeline can be run twice, with and without LCA filtering, to get two separate tracks for display.

#### Customised taxonomic filtering
More control over the taxonomic filtering can be achived using the
[scripts/rna_features/taxonomic_levels.pl](scripts/rna_features/taxonomic_levels.pl) script (see below).


#### Blacklisting and whitelisting Rfam models

Some inappropriate models may slip past the taxonomic filtering.
For the current blacklist see `rfam_blacklist` definition in [RNAFeatures_conf.pm](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/RNAFeatures_conf.pm).

Sometimes we might want to always use a model that is excluded by the default taxonomic filtering.
For the current whitelist see `rfam_whitelist` definition in [RNAFeatures_conf.pm](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/RNAFeatures_conf.pm).


#### Running the pipeline on collections (Bacteria,Fungi, Protists)

When running the pipeline on collections (ie.- bacteria, fungi, protists...)
using `-cmscan_heuristics rfam_official` reduces by a factor of 6x/7x the run-time of the CMScan analysis.
The cost in sensitivity is totally acceptable ( RFAM themselves use this parameter to run cmScan on their DB ).

The following parameters were recommended by Rfam, and they are passed automatically when choosing the  heuristic 'rfam_official':   
 * `-rfam` -- significantly increase the speed of the run without affecting much of the results quality
 * `-cut_ga option` --  automatically adapt the threshold according to the values used by the Rfam curators for each specific family, which will result in better fit cutoffs, and lead to better runtime in average. 
 * `-fit 2` -- specify the kind of output file, adding two extra columns with additional flags that help filtering the results, reducing the number of false positives.
 * `-clanin  Rfam.clanin`  -- file with groups of homologous models and some of the hits to those models will overlap (reduces false positives). Should be be downloaded from ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin 

For more details and full explanation of the Infernal software options please see:
https://docs.rfam.org/en/latest/genome-annotation.html

Example init command for bacterial collection:
``` 
CLANIN_FILE=/path/to/clanin/file # or whatever

init_pipeline.pl Bio::EnsEMBL::EGPipeline::PipeConfig::RNAFeatures_conf \
  $($CMD details script) \
  -hive_force_init 1
  -registry $REG_FILE \
  -pipeline_dir $PIPELINE_DIR \
  -eg_pipelines_dir $EG_PIPELINES_DIR  \
  -division bacteria \
  -cmscan_heuristics rfam_official \
  -production_db ensembl_production \
  -clanin_file $CLANIN_FILE \
```

#### Prepare Rfam data for release

We need to preprocess Rfam data when a new version is released, in order to format the data suitably for the pipeline.
Additionaly a file for use with customised taxonomic filters can be generated.

N.B. This should be done only once for the release.

* [scripts/rna_features/add_rfam_desc.pl](scripts/rna_features/add_rfam_desc.pl) is used to add model descritions (reported by `cmscan`) into the file with the Rfam covariance models (`Rfam.cm`). Usage instructions are given at the top of the script. Sample usage is shown below.

* [scripts/rna_features/taxonomic_levels.pl](scripts/rna_features/taxonomic_levels.pl) is used to extract data about the species that have sequences aligned for each Rfam model, cross-references with the ncbi_taxonomy database to get counts of each species (or subspecies) at the level of division and store this information into a single file. Usage instructions are given at the top of the script. Sample use cases areshown below.

Here's an example of a typical setup process:
(don't forget to get/edit/add [RFAM_VERSIONS_DIR](lib/perl/Bio/EnsEMBL/EGPipeline/PrivateConfDetails/Impl.pm.example) )

```
# initialize with proper values
RFAM_VERSION=14.2
RFAM_VERSION_DIR=... # dir to store Rfam data in

ENS_ROOT_DIR=... # dir with Ensembl repos

# get data
RFAM_DIR=${RFAM_VERSION_DIR}/${RFAM_VERSION}

mkdir -p ${RFAM_DIR}
cd ${RFAM_DIR}

lftp -e 'mirror database_files; mget README Rfam.full_region.gz Rfam.clanin Rfam.cm.gz; bye' \
  ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT

gunzip Rfam.cm.gz Rfam.full_region.gz database_files/family.txt.gz database_files/rfamseq.txt.gz

# patch Rfam with family info
perl ${ENS_ROOT_DIR}/ensembl-production-imported/scripts/rna_features/add_rfam_desc.pl \
  -cm_file Rfam.cm \
  -family_file database_files/family.txt

cmpress Rfam.cm

# generate taxonomic levels
perl ${ENS_ROOT_DIR}/ensembl-production-imported/scripts/rna_features/taxonomic_levels.pl \
  -full_region_file Rfam.full_region \
  -rfamseq_file database_files/rfamseq.txt \
  -rfam2taxonomy_file rfam2taxonomy.txt \
  > taxonomic_levels.txt

# generate database, if you need one
pushd database_files
  RFAM_VERSION_U=$(echo ${RFAM_VERSION} | perl -pe 's/\./_/') # replace '.' with '_'
  RFAM_DB="rfam_${RFAM_VERSION_U}"

  ${CMD_W} -e "drop database if exists ${RFAM_DB}; create database ${RFAM_DB};"

  ls *.sql | cut -f 1 -d . | sort | uniq |
    xargs -n 1 -I XXX \
      ${CMD_W} -D ${RFAM_DB} -e 'set foreign_key_checks = 0; source XXX.sql; set foreign_key_checks = 1'

  gunzip *txt.gz
  ls *.sql | cut -f 1 -d . | sort | uniq |
    xargs -n 1 -I XXX \
      ${CMD_W} -D ${RFAM_DB} -e 'set foreign_key_checks = 0; load data local infile "XXX.txt" into table XXX; set foreign_key_checks = 1;'

  # gzip *txt # if you need to preserve raw data
popd

# drop raw (useless) data
rm "${RFAM_DIR}/Rfam.full_region"
rm -rf "${RFAM_DIR}/database_files"  # if you don't need raw data
```

Now you can update `rfam_version` in [RNAFeatures_conf.pm](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/RNAFeatures_conf.pm).


#### Custom taxonomic filters preparation

Once you have initial `${RFAM_DIR}/taxonomic_levels.txt` and `${RFAM_DIR}/rfam2taxonomy.txt` files generated,
 you can create a custom taxonomy filter:
```
perl ${ENS_ROOT_DIR}/ensembl-production-imported/scripts/rna_features/taxonomic_levels.pl \
  -rfam2taxonomy_file ${RFAM_DIR}/rfam2taxonomy.txt \
  -levels Diptera \
  -levels Drosophila \
  -root Eukaryota \
  > /path/to/custom/taxonomic_levels/diptera.txt
```

This custom filter can be passed to init command using `-rfam_taxonomy_file` (and corresponding `-taxonomic_levels`) option, i.e.:
```
   -rfam_taxonomy_file /path/to/custom/taxonomic_levels/diptera.txt \
   -taxonomic_levels Diptera
```


### Parts
A few generic from [Common::RunnableDB](docs/Common_RunnableDB.md).

A few from [RNAFeatures](lib/perl/Bio/EnsEMBL/EGPipeline/RNAFeatures/).

