process INITIALCLASSIFICATION {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::r-base=4.2.1"                                           
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :            
        'quay.io/biocontainers/r-base:4.2.1' }"                                 

    input:
    tuple val(meta), path(hits)

    output:
    tuple val(meta), path("species.txt")   , emit: species
    tuple val(meta), path("genera.txt")    , emit: genera
    tuple val(meta), path("families.txt")  , emit: families
    tuple val(meta), env(KINGDOM)          , emit: kingdom
    tuple val(meta), env(CLASS)            , emit: classification
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    which Rscript
    Rscript --version
    Rscript --vanilla ${projectDir}/bin/sendsketch_filter.R $hits

    KINGDOM="\$(cat kingdom.txt)"
    CLASS="\$(cat classification.txt)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
