process R2TBIN {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::biopython=1.78" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'biocontainers/biopython:1.78' }"

    input:
    tuple val(meta), path(markers)

    output:
    tuple val(meta), path("${prefix}_busco_markers"), emit: busco_markers

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}_busco_markers
    r2tbinder.py ${markers} ${prefix}_busco_markers
    """
}
