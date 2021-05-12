BlastProtein

[Bio::EnsEMBL::EGPipeline::FileDump::GFF3Dumper]() is still used.
Perhaps, it's possible to replace it 
with the [Bio::EnsEMBL::Production::Pipeline::GFF3::DumpFile](https://github.com/Ensembl/ensembl-production/blob/master/modules/Bio/EnsEMBL/Production/Pipeline/GFF3/DumpFile.pm).
But `join_align_feature` param is not supported.

