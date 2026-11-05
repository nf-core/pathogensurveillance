process CALCULATE_POCP {
    tag "$group_meta.id"
    label 'process_single'

    conda "conda-forge::r-base=4.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'quay.io/biocontainers/r-base:4.2.1' }"

    input:
    tuple val(group_meta), path(gene_fam_pa)

    output:
    tuple val(group_meta), path("${prefix}_pocp.tsv"), emit: pocp
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${group_meta.id}"
    """
    calculate_pocp.R ${gene_fam_pa} ${prefix}_pocp.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
