process FILES_IN_DIR {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(dir)

    output:
    tuple val(meta), path("${dir}/*"), emit: files

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    """
}
