## EGGeneric

### Module: [Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf](../lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/EGGeneric_conf.pm).

This is a `EGPipeline`'s base class, used as a generic configuration by the other pipelines.
Contains predefined resource classes and a few options shared by other pipelines.

| option | default value |  meaning | 
| - | - | - |
`-pipeline_tag` |  | Added by some pipelines as a suffix to the eHive pipeline dbs upon creation 
`-queue_name` |  | Default LSF queue name 
`-email` | |  Default email to use

### make_resources

This method creates resource strings for LSF and SLURM, to avoid hardcoding them.

It accepts one hashref with the following keys:

| option | default value |  meaning | 
| - | - | - |
| `queue` |  | Name of the queue/partition
| `memory` |  | Memory to reserve in MB
| `time` |  | Time limit in the form 0:00:00 (h:mm:ss)
| `cpus` | 1 | number of cores
| `temp_memory` | 0 | Temp memory to reserve
| `lsf_param` | '' | A string for specific LSF parameters
| `slurm_param` | '' | A string for specific SLURM parameters
