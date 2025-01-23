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
    tuple val(group_meta), path("${prefix}_feat_seqs/cluster_*"), emit: feat_seqs
    tuple val(group_meta), path("removed_sample_ids.txt"), emit: removed_sample_ids
    tuple val(group_meta), path("removed_ref_ids.txt"), emit: removed_ref_ids


    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${group_meta.id}"
    """
    subset_busco_gene.R ${samp_ref_pairs} $min_genes $max_genes ${prefix}_feat_seqs ${busco_out_dirs}
    """
}
