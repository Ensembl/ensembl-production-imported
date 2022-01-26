## Map_interspecies_interactions
### [Bio::EnsEMBL::EGPipeline::PipeConfig::Map_interspecies_interactions_conf](../lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/Map_interspecies_interactions_conf.pm)

Storing pathogen-host interactions data from different source DBs including [PHI-base](http://www.phi-base.org) into Ensembl Interspecies interactions DB

### Prerequisites

A registry file with the urls of the interspecies_interactions_db, the ncbi_tax, and meta_db database servers specified.

A .csv file from the source_DB from whioch the interactions are sourced.

A .obo file with the descritption of the controlled terms used to describe different aspects of the interactions (optional)

### How to run


```
init_pipeline.pl Bio::EnsEMBL::EGPipeline::PipeConfig::Map_interspecies_interactions_conf \
    -pipeline_url $EHIVE_URL
    -reg_file $REGISTRY 
    -obo_file $OBO_FILE
    -source_db $SOURCE_DB
    -inputfile $INPUT_DB_FILE
```


### Parameters / Options

| option | default value |  meaning | 
| - | - | - |
| `-reg_file` |  | Different format than the one previously used in the perl eHive pipeline. This one centralises the url DB details for some databases that the pipeline is going to access. Specifically: interactions_db_url, ncbi_tax_url and meta_db_url
| `-input_file` | | path to the prepared  PHI-base current [snapshot](https://github.com/PHI-base/data/blob/master/releases/phi-base_current.csv)
| `-obo_file` | | optional path to any .obo file containing controlled vocabulary or ontologies used to describe the interactions properties
| `-source_db` | | name of the sourceDB from where we are importing data.ie.- 'PHI-base'



### Notes

N.B. The pipeline is still under construction. The first block (populating interactions to our interactionsDB) is operational but will undergo further modifications. A second part to this pipeline will extrapolate new interactions using Ensembl resources.


