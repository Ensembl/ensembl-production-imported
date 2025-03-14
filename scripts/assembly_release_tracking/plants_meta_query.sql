SELECT meta_key,meta_value FROM meta WHERE meta_key IN (
    'annotation.provider_name',
    'assembly.accession',
    'assembly.default',
    'assembly.provider_name',
    'genebuild.version',
    'ploidy',
    'species.common_name',
    'species.display_name',
    'species.scientific_name',
    'species.strain_group',
    'species.strain',
    'species.taxonomy_id',
    'strain.type'
    );