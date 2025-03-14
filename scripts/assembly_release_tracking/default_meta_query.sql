SELECT meta_key,meta_value FROM meta WHERE meta_key IN (
    'annotation.provider_name',
    'assembly.accession',
    'assembly.default',
    'assembly.provider_name',
    'genebuild.version',
    'species.common_name',
    'species.display_name',
    'species.scientific_name',
    'species.taxonomy_id'
    );