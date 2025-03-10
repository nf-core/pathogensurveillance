process MAKE_REGION_FILE {
    tag "$ref_meta.id"
    label 'process_single'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(ref_meta), path(ref)

    output:
    tuple val(ref_meta), path("*.txt"), emit: regions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${ref_meta.id}"
    """
    zgrep '>' ${ref} | sed 's/>//g' | sed 's/ .*//g' > ${prefix}.txt
    """
}
