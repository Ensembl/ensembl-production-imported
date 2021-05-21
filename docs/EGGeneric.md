## EGGeneric

### Module: [Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf](../lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/EGGeneric_conf.pm).

This is a `EGPipeline`'s base class, used as a generic configuration by the other pipelines.
Contains predefined resource classes and a few options shared by other pipelines.

| option | default value |  meaning | 
| - | - | - |
`-pipeline_tag` |  | Added by some pipelines as a suffix to the eHive pipeline dbs upon creation 
`-queue_name` |  | Default LSF queue name 
`-email` | |  Default email to use

