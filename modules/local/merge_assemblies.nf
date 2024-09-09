process MERGE_ASSEMBLIES {
    tag "All"
    label 'process_single'

    conda "conda-forge::r-base=4.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'biocontainers/r-base:4.2.1' }"

    input:
    path stats // multiple TSV outputs from FIND_ASSEMBLIES

    output:
    path "merged_assembly_stats.tsv", emit: merged_stats
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    Rscript --vanilla ${projectDir}/bin/merge_assemblies.R ${stats}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
