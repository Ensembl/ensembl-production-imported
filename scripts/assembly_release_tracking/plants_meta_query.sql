SELECT meta_key,meta_value FROM meta WHERE meta_key IN (
    'annotation.provider_name',
    'assembly.provider_name',
    'assembly.default',
    'assembly.accession',
    'genebuild.version',
    'species.scientific_name',
    'species.common_name',
    'species.taxonomy_id',
    'ploidy',
    'strain.type',
    'species.strain'
    );