process INITIAL_CLASSIFICATION {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::r-base=4.2.1"                                           
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :            
        'biocontainers/r-base:4.2.1' }"                                 

    input:
    tuple val(meta), path(hits)

    output:
    tuple val(meta), path("${prefix}_species.txt")   , emit: species
    tuple val(meta), path("${prefix}_genera.txt")    , emit: genera
    tuple val(meta), path("${prefix}_families.txt")  , emit: families
    tuple val(meta), env(KINGDOM)          , emit: kingdom
    tuple val(meta), env(CLASS)            , emit: classification
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    which Rscript
    Rscript --version
    Rscript --vanilla ${projectDir}/bin/sendsketch_filter.R $hits

    KINGDOM="\$(cat kingdom.txt)"
    CLASS="\$(cat classification.txt)"

    mv species.txt ${prefix}_species.txt 
    mv genera.txt ${prefix}_genera.txt 
    mv families.txt ${prefix}_families.txt 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
