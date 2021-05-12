# Ensembl production scripts for loading ad-hoc annotations (former eg-pipelines) 

## Prerequisites
Pipelines are intended to be run inside the Ensembl production environment.
Please, make sure you have all the proper credential, keys, etc. set up.

## Installation and configuration

### Getting this repo

```
git clone git@github.com:Ensembl/ensembl-production-imported.git
```

### Configuration

#### Refresing environment

Add `lib/perl` to `PERL5LIB` env (use instead of `modules`),
and `lib/python` to `PYTHONPATH` env 

```
export ENS_ROOT_DIR=$(pwd) # or whatever -- path to the dir to where the repo(s) was(were) cloned

export PERL5LIB=${PERL5LIB}:${ENS_ROOT_DIR}/ensembl-production-imported/lib/perl
export PYTHONPATH=${PYTHONPATH}:${ENS_ROOT_DIR}/ensembl-production-imported/lib/python
```

N.B. Please, predefine `ENS_ROOT_DIR` env.

#### Updating / setting default configuration options
To deal with the system specific configuration options 
[Bio::EnsEMBL::EGPipeline::PrivateConfDetails](lib/perl/Bio/EnsEMBL/EGPipeline/PrivateConfDetails.pm) module is used.
The actual configuration is loaded from [Bio::EnsEMBL::EGPipeline::PrivateConfDetails::Impl](lib/perl/Bio/EnsEMBL/EGPipeline/PrivateConfDetails/Impl.pm.example).

All the used options are listed in [Impl.pm.example](lib/perl/Bio/EnsEMBL/EGPipeline/PrivateConfDetails/Impl.pm.example).
Please, define them before running pipelines.

This can be done either by copying this file and editing it.
```
cp ensembl-production-imported/lib/perl/Bio/EnsEMBL/EGPipeline/PrivateConfDetails/Impl.pm{.example,}
# edit ensembl-production-imported/lib/perl/Bio/EnsEMBL/EGPipeline/PrivateConfDetails/Impl.pm
```

Or by creating a separate repo with `lib/perl/Bio/EnsEMBL/EGPipeline/PrivateConfDetails/Impl.pm` and adding corresponding `lib/perl` to your `PERL5LIB` env.

#### Resources and queues

You can override the default queue used to run pipeline by adding
`-queue_name` option to the `init_pipeline.pl` command (see below).

### Initialising and running pipelines


Every pipeline is derived from 
[Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/EGGeneric_conf.pm)
(see [EGGeneric documentation](docs/EGGeneric.md)) for details.

And the same perl class prefix used for every pipeline:
  `Bio::EnsEMBL::EGPipeline::PipeConfig::` .

N.B. Don't forget to specify `-reg_file` option for the `beekeeper.pl -url $url -reg_file $REG_FILE -loop` command.

```
init_pipeline.pl Bio::EnsEMBL::EGPipeline::PipeConfig::RNAFeatures_conf \
    $($CMD details script) \
    -hive_force_init 1\
    -queue_name $SPECIFIC_QUEUE_NAME \
    -registry $REG_FILE \
    -production_db "$($PROD_SERVER details url)""$PROD_DBNAME" \
    -pipeline_tag "_${SPECIES_TAG}" \
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


## Pipelines
| Pipeline name | Module | Description | Document | Comment|
| - | - | - | - | - |
| EGGeneric | [Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/EGGeneric_conf.pm) | generic pipeline configuration  | [EGGeneric](docs/EGGeneric.md) | All other pipelines are derived from this one
| RepeatModeler | [Bio::EnsEMBL::EGPipeline::PipeConfig::RepeatModeler_conf](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/RepeatModeler_conf.pm) | Building de-nove repeat libs | [RepeatModeler](docs/RepeatModeler.md) | | 
| DNAFeatures | [Bio::EnsEMBL::EGPipeline::PipeConfig::DNAFeatures_conf](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/DNAFeatures_conf.pm) | repeat masking  | [DNAFeatures](docs/DNAFeatures.md) | redat_repeatmasker_library  should be explicitly specified |
| RNAFeatures | [Bio::EnsEMBL::EGPipeline::PipeConfig::RNAFeatures_conf](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/RNAFeatures_conf.pm) | Non-coding rna features (tRNA, miRNA, etc) discovery  | [RNAFeatures](docs/RNAFeatures.md) |
| RNAGenes | [Bio::EnsEMBL::EGPipeline::PipeConfig::RNAGenes_conf](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/RNAGenes_conf.pm) | Create non-coding genes based on rna features | [RNAGenes](docs/RNAGenes.md) | Specify id_db_{host,port,user,dbname,...} options if run_context != "VB"
| SRAAlignment_BRC4 | [Bio::EnsEMBL::EGPipeline::PipeConfig::SRAAlignment_BRC4_conf](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/SRAAlignment_BRC4_conf.pm) | Perform RNA(DNA) short read aligments | [SRAAlignment_BRC4](docs/SRAAlignment_BRC4.md) |
| WGA2GenesDirect | [Bio::EnsEMBL::EGPipeline::PipeConfig::WGA2GenesDirect_conf](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/WGA2GenesDirect_conf.pm) | Project transripts and create genes based on compara lastz mappings | [WGA2GenesDirect](docs/WGA2GenesDirect.md)
| Xref_GPR | [Bio::EnsEMBL::EGPipeline::PipeConfig::Xref_GPR_conf](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/Xref_GPR_conf.pm) | Load Plant Reactome data | [Xref_GPR](docs/Xref_GPR.md) | use `-uppercase_gene_id 1` option to allow usage of uppercase gene stable IDs for mapping (i.e. for _Oryza sativa_ (rice))
| AlignmentXref | [Bio::EnsEMBL::EGPipeline::PipeConfig::AlignmentXref_conf](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/AlignmentXref_conf.pm) | Alignment bases xrefs | [AlignmentXref](docs/AlignmentXref.md) | Used as a part of the `AllXref` pipeline
| Xref | [Bio::EnsEMBL::EGPipeline::PipeConfig::Xref_conf](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/Xref_conf.pm) | MD5-based UniParc/Uniprot Xref pipeline  | [Xref](docs/Xref.md) | Used as a part of the `AllXref` pipeline
| AllXref | [Bio::EnsEMBL::EGPipeline::PipeConfig::AllXref_conf](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/AllXref_conf.pm) | Combined Xref/AlignmentXref pipeline | [AllXref](docs/AllXref.md) |
| FindPHIBaseCandidates | [Bio::EnsEMBL::EGPipeline::PipeConfig::FindPHIBaseCandidates_conf](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/FindPHIBaseCandidates_conf.pm) | Load Xrefs from PHIBase | [FindPHIBaseCandidates](docs/FindPHIBaseCandidates.md)

## Obsolete pipelines
| Pipeline name | Module | Description | Document | Comment | Alternative |
| - | - | - | - | - | - |
| AnalyzeTables | [Bio::EnsEMBL::EGPipeline::PipeConfig::AnalyzeTables_conf](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/AnalyzeTables_conf.pm) | Runs SQL ANALIZE / OPTIMIZE on tables for DBs present in the registry
| EC2Rhea | [Bio::EnsEMBL::EGPipeline::PipeConfig::EC2Rhea_conf](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/EC2Rhea_conf.pm) | Adding chemical and transport reactions (Rhea2RC) xrefs (used by 'microbes')  | | Specify `ec2rhea_file` as there's no default |
| ExonerateAlignment | [Bio::EnsEMBL::EGPipeline::PipeConfig::ExonerateAlignment_conf](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/ExonerateAlignment_conf.pm) | Aligning Fasta files to a genome with Exonerate  |  | Specify `-exonerate_2_4_dir` option if use _exonerate-server_ (` -use_exonerate_server 1`)
| ShortReadAlignment | [Bio::EnsEMBL::EGPipeline::PipeConfig::ShortReadAlignment_conf](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/ShortReadAlignment_conf.pm)
| STARAlignment | [Bio::EnsEMBL::EGPipeline::PipeConfig::STARAlignment_conf](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/STARAlignment_conf.pm) 
| BlastNucleotide | [Bio::EnsEMBL::EGPipeline::PipeConfig::BlastNucleotide_conf](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/BlastNucleotide_conf.pm)
| BlastProtein | [Bio::EnsEMBL::EGPipeline::PipeConfig::BlastProtein_conf](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/BlastProtein_conf.pm) | | | `EGPipeline::FileDump::GFF3Dumper` could not be replaced with [Production::Pipeline::GFF3::DumpFile](https://github.com/Ensembl/ensembl-production/blob/master/modules/Bio/EnsEMBL/Production/Pipeline/GFF3/DumpFile.pm) as no `join_align_feature` param is provided | | 
| Bam2BigWig | [Bio::EnsEMBL::EGPipeline::PipeConfig::Bam2BigWig_conf](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/Bam2BigWig_conf.pm) |
| ProjectGenes | [Bio::EnsEMBL::EGPipeline::PipeConfig::ProjectGenes_conf](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/ProjectGenes_conf.pm)
| ProjectGeneDesc | [Bio::EnsEMBL::EGPipeline::PipeConfig::ProjectGeneDesc_conf](lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/ProjectGeneDesc_conf.pm)

## Replaced pipelines
| Old pipeline module | Alternative | Description | Document | Comment |
| - | - | - | - | - |
| CoreStatistics | [Bio::EnsEMBL::Production::Pipeline::PipeConfig::CoreStatistics_conf](https://github.com/Ensembl/ensembl-production/blob/master/modules/Bio/EnsEMBL/Production/Pipeline/PipeConfig/CoreStatistics_conf.pm) | Core stats pipeline | | use `-skip_metadata_check 1` if core is not submitted (always for new species) 
| FileDump | [Bio::EnsEMBL::Production::Pipeline::PipeConfig::FileDump_conf](https://github.com/Ensembl/ensembl-production/blob/master/modules/Bio/EnsEMBL/Production/Pipeline/PipeConfig/FileDump_conf.pm) | Serialize core 
| FileDump{Compara,GFF} | same as above
| FileDumpVEP | [Bio::EnsEMBL::Production::Pipeline::PipeConfig::FileDumpVEP_conf](https://github.com/Ensembl/ensembl-production/blob/master/modules/Bio/EnsEMBL/Production/Pipeline/PipeConfig/FileDumpVEP_conf.pm) | Dump VEP data
| LoadGFF3 | [Bio::EnsEMBL::Pipeline::PipeConfig::LoadGFF3_conf](https://github.com/MatBarba/new_genome_loader/blob/master/lib/perl/Bio/EnsEMBL/Pipeline/PipeConfig/LoadGFF3_conf.pm) | Load gene models from GFF3 and accompanied files | See [new_genome_loader](https://github.com/MatBarba/new_genome_loader) for details
| LoadGFF3Batch | [Bio::EnsEMBL::Pipeline::PipeConfig::LoadGFF3Batch_conf](https://github.com/MatBarba/new_genome_loader/blob/master/lib/perl/Bio/EnsEMBL/Pipeline/PipeConfig/LoadGFF3Batch_conf.pm) | Batch load models from GFF3 files  | See [new_genome_loader](https://github.com/MatBarba/new_genome_loader) for details
| GeneTreeHighlighting | [Bio::EnsEMBL::Production::Pipeline::PipeConfig::GeneTreeHighlighting](https://github.com/Ensembl/ensembl-production/blob/master/modules/Bio/EnsEMBL/Production/Pipeline/PipeConfig/GeneTreeHighlighting_conf.pm)  | Populate compara table with GO and InterPro terms, to enable highlighting |
| GetOrthologs | [Bio::EnsEMBL::Production::Pipeline::PipeConfig::DumpOrtholog](https://github.com/Ensembl/ensembl-production/blob/master/modules/Bio/EnsEMBL/Production/Pipeline/PipeConfig/DumpOrtholog_conf.pm)



## Runnables worth additional mentioning
| Runnable |  Description | Document | Comment|
| - | - | - | - |
| Common::RunnableDB::CreateOFDatabase | 
| Analysis::Config::General | 

## Scripts

| Script |  Description | Document | Comment|
| - | - | - | - |
| [brc4/repeat_for_masker.pl](scripts/brc4/repeat_for_masker.pl) | ....
| [brc4/repeat_tab_to_list.pl](scripts/brc4/repeat_tab_to_list.pl) | ....
| [misc_scripts/get_trans.pl](scripts/misc_scripts/get_trans.pl) | get transcriptions and tranaslations | In pipelines use Bio::EnsEMBL::EGPipeline::Common::RunnableDB::DumpProteome  and Bio::EnsEMBL::EGPipeline::Common::RunnableDB::DumpTranscriptome 
| [misc_scripts/load_xref.pl](scripts/misc_scripts/load_xref.pl)
| [misc_scripts/remove_entities.pl](scripts/misc_scripts/remove_entities.pl)
| [misc_scripts/gene_stable_id_mapping.pl](scripts/misc_scripts/gene_stable_id_mapping.pl)
| [misc_scripts/add_karyotype.pl](scripts/misc_scripts/add_karyotype.pl)
| [misc_scripts/load_karyotype_from_gff.pl](scripts/misc_scripts/load_karyotype_from_gff.pl) 
| [misc_scripts/gene_stable_id_mapping.pl](scripts/misc_scripts/gene_stable_id_mapping.pl) 
| [rna_features/add_rfam_desc.pl](scripts/rna_features/add_rfam_desc.pl) | prepare Rfam db for RNAFeatures | [RNAFeatures](docs/RNAFeatures.md)
| [rna_features/taxonomic_levels.pl](scripts/rna_features/taxonomic_levels.pl) | prepare Rfam db for RNAFeatures | [RNAFeatures](docs/RNAFeatures.md)
| [phi_ontology/phi-base_ontologies.pl](scripts/phi_ontology/phi-base_ontologies.pl) | normalising phi-base data .csv based on onlologies in [scripts/phi_ontology](scripts/phi_ontology) | [FindPHIBaseCandidates](docs/FindPHIBaseCandidates.md)

## Replaced scripts

| Script |  Substitution | Document | Comment|
| - | - | - | - |
| production_db/analysis_desc_from_prod.pl | [Bio::EnsEMBL::Production::Pipeline::PipeConfig::ProductionDBSync_conf](https://github.com/Ensembl/ensembl-production/blob/master/modules/Bio/EnsEMBL/Production/Pipeline/PipeConfig/ProductionDBSync_conf.pm)
| production_db/attrib_type_from_prod.pl | [Bio::EnsEMBL::Production::Pipeline::PipeConfig::ProductionDBSync_conf](https://github.com/Ensembl/ensembl-production/blob/master/modules/Bio/EnsEMBL/Production/Pipeline/PipeConfig/ProductionDBSync_conf.pm)
| production_db/external_db_from_prod.pl | [Bio::EnsEMBL::Production::Pipeline::PipeConfig::ProductionDBSync_conf](https://github.com/Ensembl/ensembl-production/blob/master/modules/Bio/EnsEMBL/Production/Pipeline/PipeConfig/ProductionDBSync_conf.pm)
| production_db/add_species_analysis.pl | [Bio::EnsEMBL::Production::Pipeline::PipeConfig::ProductionDBSync_conf](https://github.com/Ensembl/ensembl-production/blob/master/modules/Bio/EnsEMBL/Production/Pipeline/PipeConfig/ProductionDBSync_conf.pm)


## Various docs
See [docs](docs)

## TODO
Tests, tests, tests...

## Acknowledgements

For obvoius reason the whole history of the source project had to go.
Most of this code and documentation is inherited from the [EnsemblGenomes](https://github.com/EnsemblGenomes) project.

We appreciate the effort and time spent by developers of the [EnsemblGenomes](https://github.com/EnsemblGenomes) project.

Thank you!

