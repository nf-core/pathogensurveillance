process FILTER_ASSEMBLY {
    tag "$meta.id"

    conda (params.enable_conda ? "conda-forge::biopython=1.78" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'biocontainers/biopython:1.78' }"

    input:
    tuple val(meta), path(scaffold)

    output:
    tuple val(meta), path("*_filtered.fasta"), emit: filtered
    path "versions.yml"             , emit: versions

    script:
    def args = task.ext.args ?: ''                                              
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gzip -dc ${scaffold} > ${prefix}_unzipped.fasta

    process_spades_assem.py \\
        --verbose \\
        --summary ${prefix}.summary \\
        $args \\
        ${prefix}_unzipped.fasta > ${prefix}_filtered.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
