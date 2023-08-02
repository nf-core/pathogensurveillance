process ASSIGN_GROUP_REFERENCES {
    tag "$group_meta.id"
    label 'process_low'

    conda "conda-forge::r-base=4.2.1"                                           
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :            
        'quay.io/biocontainers/r-base:4.2.1' }"                                 

    input:
    tuple val(group_meta), path(ani_matrix), path(samp_ref_pairs)

    output:
    tuple val(group_meta), path("${prefix}.csv"), emit: samp_ref_pairs
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${group_meta.id}"
    """
    assign_group_reference.R ${ani_matrix} ${samp_ref_pairs} ${prefix}.csv

    cat <<-END_VERSIONS > versions.yml                                          
    "${task.process}":                                                          
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS                                                                
    """
}
