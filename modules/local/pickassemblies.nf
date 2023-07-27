process PICK_ASSEMBLIES {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::r-base=4.2.1"                                           
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :            
        'quay.io/biocontainers/r-base:4.2.1' }"                                 

    input:
    tuple val(meta), path(families), path(genera), path(species)
    path stats

    output:
    tuple val(meta), path("${prefix}.tsv")     , emit: stats
    tuple val(meta), path("${prefix}_ids.txt"), emit: id_list
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    pick_assemblies.R ${families} ${genera} ${species} ${stats} 5 ${prefix}.tsv

    tail -n +2 ${prefix}.tsv | cut -f2 > ${prefix}_ids.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
