process EXTRACT_PIRATE_GENES {
    tag "$group_meta.id"
    label 'process_low'

    conda "conda-forge::r-base=4.2.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'biocontainers/r-base:4.2.1' }"

    input:
    tuple val(group_meta), path(pirate_results, stageAs: 'pirate_results/*')

    output:
    tuple val(group_meta), path("gene_sequences"), emit: gene_sequences
    path "versions.yml"                            , emit: versions_extract_pirate_genes, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${group_meta.id}"
    """
    extract_pirate_genes.R pirate_results .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
