process SUBSET_CORE_GENES {
    tag "$group_meta.id"
    label 'process_single'

    conda "conda-forge::r-base=4.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'biocontainers/r-base:4.2.1' }"

    input:
    tuple val(group_meta), path(gene_fam), path(feat_seqs), path(sample_data)
    val min_core_genes
    val max_core_genes

    output:
    tuple val(group_meta), path("core_genes/${prefix}--cluster_*.tsv"), emit: gene_fam
    tuple val(group_meta), path("feat_seqs/${prefix}--cluster_*")     , emit: feat_seq
    tuple val(group_meta), path("message_data.tsv")                              , emit: message_data
    path "versions.yml"                                                          , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${group_meta.id}"
    """
    subset_core_gene.R $gene_fam $feat_seqs $sample_data $min_core_genes $max_core_genes ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
