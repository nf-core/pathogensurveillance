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
    val n_ref_strains
    val n_ref_species
    val n_ref_genera

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: stats
    path "merged_assembly_stats.tsv"      , emit: merged_stats
    tuple val(meta), env(COUNT)           , emit: line_count
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    pick_assemblies.R ${families} ${genera} ${species} ${n_ref_strains} ${n_ref_species} ${n_ref_genera} ${prefix}.tsv ${assem_data_tsvs}
    COUNT=\$(cat ${prefix}.tsv | wc -l)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
