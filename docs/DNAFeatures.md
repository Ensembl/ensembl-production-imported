## DNAFeatures
### Module: [Bio::EnsEMBL::EGPipeline::PipeConfig::DNAFeatures_conf](../lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/DNAFeatures_conf.pm)

The DNA Features pipeline runs three programs to annotate repeat features: RepeatMasker (repeatmask), DustMasker (dust), and TRF (trf).

### Prerequisites
A registry file with the locations of the core database server(s) and the production database (or `-production_db $PROD_DB_URL` specified).

### How to run
```
# REP_LIB_OPT= # or whatever

init_pipeline.pl Bio::EnsEMBL::EGPipeline::PipeConfig::DNAFeatures_conf \
  $($CMD details script) \
  -registry $REG_FILE \
  -production_db "$($PROD_SERVER details url)""$PROD_DBNAME" \
  -hive_force_init 1\
  -pipeline_tag "_${SPECIES_TAG}" \
  -pipeline_dir $OUT_DIR \
  -species $SPECIES \
  -redatrepeatmasker 0 \
  -always_use_repbase 1 \
  -repeatmasker_timer '10H' \
  $REP_LIB_OPT \
  -repeatmasker_repbase_species "$REPBASE_SPECIES_NAME" \
  -max_seq_length 300000 \
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
| `-division` | | division (intersection with registry) to be processed 
| `-species` |  | species to process, several `-species ` options are possible
| `-pipeline_dir` | | directory to store results to
| `-report_dir ` | | direcory to write reports to
| `-dust` | 1 | run DustMasker (dust); 0 -- not to run 
| `-trf` | 1 | run TRF(trf); 0 -- not to run
| `-repeatmasker` | 1 | run RepeatMasker (repeatmask); 0 -- not to run.
| `-redatrepeatmasker` | 0 | use plant-specific library of DNA repeats [REdat](https://pgsb.helmholtz-muenchen.de/plant/recat/); 1 -- to run
| `-redat_repeatmasker_library` | [REDAT_REPEATMASKER_LIBRARY_PATH](../lib/perl/Bio/EnsEMBL/EGPipeline/PrivateConfDetails/Impl.pm.example) | path to REdat library, to be used with `-redatrepeatmasker 1`
| `-repeatmasker_library` | custom RepeatMasker libraries; several are allowed, of the form `-repeatmasker_library ${SPECIES_NAME_1}=/path/to/lib/file`
| `-repeatmasker_logic_name` | | use custom (pre-existing!) logic name instead of `repeatmask_customlib` when `-repeatmasker_library` is used (see notes below).
| `-always_use_repbase` | 0 | 1 -- to always run default (RepBase one) `repeatmasker`, even if custom `-repeatmasker_library` specifed (see notes below)
| `-repeatmasker_repbase_species` | same as `-species` | use this name to run default (RepBase-bases) repeatmasker   
| `-parameters_json` | | JSON with species specific parameters (see below)
| `-repeatmasker_sensitivity ${SPECIES}\|_all_=<level>` | `all=automatic` | controls RepeatMasker trade-off between sensitivity and run-time. Possible levels: `automatic`, `very_low`, `low`, `medium`, `high`, `very_high`. See below.
| `-delete_existing` | 1 | Delete preexisting `repeatmasker` analysis data; 0 to disable
| `-max_seq_length` | 1_000_000 | maximum length for a single scaffold, scaffolds will be splitted on exciding this value
| `-max_hive_capacity` | 50 | The default hive capacity (i.e. the maximum number of concurrent workers) (no enforced upper limit, but above 150 you might run into problems with locking or database connections).
| `-repeatmasker_timer` | `16H` | default time limit for the RepeatMasker analysis 
| `-pipeline_name` | `dna_features_${ENS_VERSION}_${pipeline_tag}` | The hive database name will be `${USER}_${pipeline_name}`
| `-pipeline_tag` |  | Tag to append to the  default `-pipeline_name`
| `-repeatmasker_exe` | RepeatMasker | Path to the RepeatMasker executable file.
| `-repeatmasker_parameters` | `-nolow -s -gccalc` | Default is to exclude low-complexity annotations; use slower (and more sensitive) search; and calculate (rather than estimate) GC content 
| `-dust_exe` | dustmasker | Path to the dust executable file.
| `-dust_parameters_hash` | `{}` | DustMasker module takes a hash of parameters, rather than requiring explicit command line options like RepeatMasker. It is generally not necessary to override the default parameters.
| `-trf_exe` | trf | Path to the TRF executable file.
| `-trf_parameters_hash` | `{}` | TRF module (in ensembl-analysis) takes a hash of parameters, rather than requiring explicit command line options like RepeatMasker. It is generally not necessary to override the default parameter. 
| `-production_lookup` | 1 |  Fetch analysis display name, description and web data from the production database; 0 -- to disable
| `-email_repeat_report` | 1 | Send an email with a summary of the repeats when the pipeline finishes; 0 -- to disable


#### Notes

To add use plants specific *nrTEplants* libary, define `$NRPLANTSLIB` path and add this ooption to the init command:
```
  -repeatmasker_library zea_mays=$NRPLANTSLIB \
  -repeatmasker_logic_name zea_mays=repeatmask_nrplants
```

To run both custom repeatmask and default (RepBase) one 
```
  -repeatmasker_library bombyx_mori=/path/to/lib/file \
  -repeatmasker_library pediculus_humanus=/path/to/another/lib/file \
  -always_use_repbase 1
```

To use different species name for default (RepBase-based) repeatmasker run
```
  -species dinothrombium_tinctorium \
  -repeatmasker_repbase_species Tetranychus_urticae \
  -repeatmasker_library dinothrombium_tinctorium=/path/to/another/lib/file \
  -always_use_repbase 1
```

#### Species specific parameters JSON

If you need to run the pipeline for a lot of species,
it might be preferable to load the species specific parameters from a file.
To do so, create a json file as a dictionary with the keys being the parameter names, like the following example.
```
{
   "species" : [
      "pediculus_humanus",
   ]
   "repeatmasker_library" : {
      "pediculus_humanus" : "/path/to/lib/pediculus_humanus.lib",
   },
}
```
.

And load it with `--parameters_json <path/to/json>` when running init command.

#### Controlling RepeatMasker sensitivity and run-time

There's a trade-off with RepeatMasker between sensitivity and run-time.
Sensitivity is controlled via the alignment engine (crossmatch or a variant of NCBI BLAST+) and some command line flags,
but it's a pain for users of this pipeline to have to think about these sorts of things;
so, left to its own defaults, the pipeline will do its best to use sensible settings.
This default behaviour can be overridden, on either a species-specific or pipeline-wide level,
via the `-repeatmasker_sensitivity` parameter, which can be set to one of:
 `automatic` (default), `very_low`, `low`, `medium`, `high`, `very_high`.
The rationale and workings for the automatic settings are detailed elsewhere, but briefly,
it takes into account the genome size, the size of each DNA chunk, and the size of the repeat library, and sets sensitivity at `low`, `medium` or `high`.

The `very_high` option is only recommended if you absopositively need to wring every ounce of sensitivity out of RepeatMasker; compared to `high`, the coverage gains will be marginal (fractions of a percentage point), and the execution time will be exponentially longer. Conversely, there's only a tiny amount of time to be gained by using `very_low` as opposed to `low`, and there's not likely to be much difference (if any) in sensitivity, so the `very_low` option is probably only useful for running against a massive number of genomes (e.g. bacteria).

If running for a custom library with the `-always_use_repbase` parameter,
any specificity parameter will apply to both analyses;
the only way to avoid this is to run the pipeline twice, once for RepBase, once for the custom library.

An example set of options to be added to run one species with high sensitivity and let the other be automatically determined:
```
  -species pediculus_humanus \
  -repeatmasker_library pediculus_humanus=/path/to/another/lib/file \
  -repeatmasker_sensitivity pediculus_humanus=high
```

To force a particular sensitivity level across all species, use 'all' instead of the species name:
```
  -repeatmasker_sensitivity all=low
```

### Parts
A few generic from [Common::RunnableDB](../docs/Common_RunnableDB.md).

A few from [DNAFeatures](../lib/perl/Bio/EnsEMBL/EGPipeline/DNAFeatures/).

