process ASSIGN_CONTEXT_REFERENCES {
    tag "$group_meta.id"
    label 'process_single'

    conda "conda-forge::r-base=4.2.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'biocontainers/r-base:4.2.1' }"

    input:
    tuple val(group_meta), path(ani_matrix), val(samp_ref_pairs_raw)
    val n_ref_closest
    val n_ref_closest_named
    val n_ref_context

    output:
    tuple val(group_meta), path("${prefix}_context_refs.tsv"), emit: references
    path "versions.yml"                                      , emit: versions_assign_context_references, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${group_meta.id}"
    """
    echo "${samp_ref_pairs_raw}" | base64 -d > ${prefix}_samp_ref_pairs.tsv

    assign_context_references.R ${ani_matrix} ${prefix}_samp_ref_pairs.tsv ${n_ref_closest} ${n_ref_closest_named} ${n_ref_context} ${prefix}_context_refs.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
