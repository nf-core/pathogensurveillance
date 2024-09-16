process R2TDIR {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(busco_dir)

    output:
    tuple val(meta), path("${prefix}_r2t_markers"), emit: markers
    tuple val(meta), path("${prefix}_dna_ref.fa"), emit: dna_ref

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}_r2t_markers
    cat **/*ntformatted.fa > ${prefix}_dna_ref.fa
    cp **/*aaformatted.fa ${prefix}_r2t_markers/
    """
}
