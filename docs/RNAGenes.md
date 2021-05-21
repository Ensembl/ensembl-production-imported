## RNAGenes
#### Module [Bio::EnsEMBL::EGPipeline::PipeConfig::RNAGenes_conf](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/RNAGenes_conf.pm)

This pipeline annotates RNA genes based on Rfam alignments, tRNAscan predictions, and miRBase data.
A subset of the alignments produced by the [RNAFeatures_conf](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/RNAFeatures_conf.pm) ([doc](docs/RNAFeatures.md)),
with strict taxonomic filtering, is used as a requisite.


### Prerequisites

A registry file with the locations of the core database server(s) and the production database (or `-production_db $PROD_DB_URL` specified).


### How to run

```
init_pipeline.pl Bio::EnsEMBL::EGPipeline::PipeConfig::RNAGenes_conf \
    $($CMD details script) \
    -hive_force_init 1\
    -registry $REG_FILE \
    -production_db "$($PROD_SERVER details url)""$PROD_DBNAME" \
    -pipeline_tag "_${SPECIES_TAG}" \
    -pipeline_dir $OUT_DIR/rna_genes \
    -species $SPECIES \
    -eg_pipelines_dir $ENS_DIR/ensembl-production-imported \
    -all_new_species 1 \
    -run_context vb \
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
| `-registry` | | registry file with the locations of the core DBs and production DB
| `-production_db` | | connection URL for the production DB; provide if no production DB is in the registry
| `-old_registry` | | registry file with the locations of the core DBs for the previous release; no need if annotating only new genomes
| `-all_new_species` | 0 |  Stop doing stable ID, all species are assumed to be new; 1 -- to disable mapping
| `-use_cmscan` | 1 | create genes from `cmscan` (Rfam mostly) alignments; 0 -- to skip
| `-use_trnascan` | 1 | create genes from `tRNAscan` alignments; 0 -- to skip
| `-use_mirbase` | 1 | create genes from `mirBase` alignments; 0 -- to skip
| `-run_context` | `eg` | style of the stable identifiers: `eg` -- `ENSRNA\d{9}` like; `vb` -- use `species.stable_id_prefix` from the core DB metatable 
| `-gene_source` | `species.division` from meta if defined; `Ensembl` -- otherwise | name to use as a gene source
| `-mirbase_source_logic_name` |  `mirbase` | `logic_name` of the source alignments (already existing in the DB) 
| `-mirbase_target_logic_name` | `mirbase_gene` | `logic_name` for the genes to be created
| `-trnascan_source_logic_name` | `trnascan_align` | `logic_name` of the source alignments (already existing in the DB)
| `-trnascan_target_logic_name` | `trnascan_gene` | `logic_name` for the genes to be created
| `-rfam_version` | [RFAM_VERSION](lib/perl/Bio/EnsEMBL/EGPipeline/PrivateConfDetails/Impl.pm.example) | set `rfam_version` to be included into default `-cmscan_source_logic_name` and `-cmscan_target_logic_name` values
| `-cmscan_source_logic_name` | `cmscan_rfam_${rfam_version}_lca` | `logic_name` of the source alignments (default `_lca` assumes strict taxonomic filtering was used)
| `-cmscan_target_logic_name` | `rfam_${rfam_version}_gene` | `logic_name` for the genes to be created
| `-id_db_host` | [ENSEMBL_ENA_IDENTIFIERS_HOST](lib/perl/Bio/EnsEMBL/EGPipeline/PrivateConfDetails/Impl.pm.example) | connection details for the stable IDs DB
| `-id_db_port` | [ENSEMBL_ENA_IDENTIFIERS_PORT](lib/perl/Bio/EnsEMBL/EGPipeline/PrivateConfDetails/Impl.pm.example) | connection details for the stable IDs DB
| `-id_db_user` | [ENSEMBL_ENA_IDENTIFIERS_USER](lib/perl/Bio/EnsEMBL/EGPipeline/PrivateConfDetails/Impl.pm.example) | connection details for the stable IDs DB
| `-id_db_dbname` | [ENSEMBL_ENA_IDENTIFIERS_DBNAME](lib/perl/Bio/EnsEMBL/EGPipeline/PrivateConfDetails/Impl.pm.example) | connection details for the stable IDs DB
| `-id_db_pass` | | connection details for the stable IDs DB
| `-pipeline_tag` |  | Tag to append to the  default `-pipeline_name`
| `-pipeline_name` | `rna_genes_${ENS_VERSION}_<pipeline_tag>` | The hive database name will be `${USER}_${pipeline_name}`
| `-production_lookup` | 1 |  Fetch analysis display name, description and web data from the production database; 0 -- to disable


### Notes

#### Filtering
The taxonomic filtering itself, if enabled, for the [RNAFeatures_conf](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/RNAFeatures_conf.pm) ([doc](docs/RNAFeatures.md)), removes many false positives.
This pipeline assumes that RNA gene creation is conservative,
and expects the taxonomic filtering be based on shared ancestry (rather than at the level of divisions).

A few filters applied automatically:
 * regulatory elements (everything with "misc_RNA" biotype) are filtered out;
 * palindromic RNA sequences can lead to overalpping credible hits (on the same or different strand),
  only one such alignment is picked to be converted to a gene, based on the better E-value.

`miRBase` alignments are trusted and not filtered.

`cmscan` filters can be tweaked using these pipeline options (see [Inferal docs](http://eddylab.org/infernal/Userguide.pdf) for some details):

| option | default value |  meaning | 
| - | - | - |
| `-evalue_threshold` | 1e-6 | `cmscan` E-value threshold
| `-truncated` | 0 | allow processing of truncated alignments; 1 -- partial genes are allowed
| `-nonsignificant` | 0 | allow processing of  non-significant alignments; 1 -- to allow
| `-bias_threshold` | 0.3 | maximum degree of allowable GC/AT bias
| `-has_structure` | 1 | only use features if they have structure; 0 -- allow to use features lacking structure 
| `-allow_repeat_overlap` | 1 | allow genes which overlap a repeat feature; 0 -- to disallow (not recommended)
| `-allow_coding_overlap` | 0 | allow genes which overlap a protein-coding exon; 1 -- to allow
| `-maximum_per_hit_name` | `{ pre_miRNA' => 100, }` | limits on usage/sharing of the same hit(model) name

`tRNAscan` filters: 
 * always applied:
   * tRNA genes are not allowed to overlap repeat regions
   * tRNA genes are not allowed to overlap protein-coding exons

Configurable:

| option | default value |  meaning | 
| - | - | - |
| `-score_threshold` |  40 | a threshold for the COVE score


#### Preserving stable IDs between re-runs / releases

There's an option to preserve stable IDs between releases and reruns runs using `-id_db_*` database options. 


### Parts
A few generic from [Common::RunnableDB](docs/Common_RunnableDB.md).

A few from [RNAFeatures](lib/perl/Bio/EnsEMBL/EGPipeline/RNAFeatures/).

