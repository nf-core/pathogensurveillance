process SUBSET_SHARED_GENES {
    tag "$group_meta.id"
    label 'process_single'

    conda "conda-forge::r-base=4.2.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'biocontainers/r-base:4.2.1' }"

    input:
    tuple val(group_meta), path(gene_sequences), val(metadata)
    val min_genes
    val max_genes

    output:
    tuple val(group_meta), path("feat_seqs/${prefix}--cluster_*")     , emit: feat_seqs, optional: true
    tuple val(group_meta), path("${prefix}--cluster_*_gene_list.tsv") , emit: gene_lists, optional: true
    tuple val(group_meta), path("removed_ids.tsv")                    , emit: removed_ids, optional: true
    path "versions.yml"                                               , emit: versions_subset_shared_genes, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${group_meta.id}--${metadata.source}"
    """
    echo "${metadata.samp_ref_pairs}" | base64 -d > ${group_meta.id}_samp_ref_pairs.tsv

    subset_shared_genes.R ${gene_sequences} ${group_meta.id}_samp_ref_pairs.tsv ${min_genes} ${max_genes} ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
