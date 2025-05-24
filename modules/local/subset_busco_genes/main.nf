process SUBSET_BUSCO_GENES {
    tag "$group_meta.id"
    label 'process_single'

    conda "conda-forge::r-base=4.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'biocontainers/r-base:4.2.1' }"

    input:
    tuple val(group_meta), path(busco_out_dirs), path(samp_ref_pairs)
    val min_genes
    val max_genes

    output:
    tuple val(group_meta), path("${prefix}_feat_seqs/${prefix}--cluster_*"), emit: feat_seqs
    tuple val(group_meta), path("message_data.tsv")             , emit: message_data
    path "versions.yml"                                         , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${group_meta.id}"
    """
    subset_busco_gene.R ${samp_ref_pairs} $min_genes $max_genes ${prefix} ${busco_out_dirs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
