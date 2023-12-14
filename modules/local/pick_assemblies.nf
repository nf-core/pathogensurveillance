process PICK_ASSEMBLIES {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::r-base=4.2.1"                                           
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :            
        'quay.io/biocontainers/r-base:4.2.1' }"                                 

    input:
    tuple val(meta), path(families), path(genera), path(species)
    path assem_data_tsvs

    output:
    tuple val(meta), path("${prefix}.tsv")    , emit: stats
    tuple val(meta), path("${prefix}_ids.txt"), emit: id_list
    path "merged_assembly_stats.tsv", emit: merged_stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    Rscript --vanilla ${projectDir}/bin/pick_assemblies.R ${families} ${genera} ${species} 5 ${prefix}.tsv ${assem_data_tsvs}

    tail -n +2 ${prefix}.tsv | cut -f1,3 > ${prefix}_ids.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
