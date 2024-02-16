process SUBSET_CORE_GENES {
    tag "$ref_meta.id"
    label 'process_single'

    conda "conda-forge::r-base=4.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'biocontainers/r-base:4.2.1' }"

    input:
    tuple val(ref_meta), path(gene_fam), path(feat_seqs)
    path sample_data
    val min_core_genes
    val min_core_samps
    val min_core_refs
    val max_core_genes

    output:
    tuple val(ref_meta), path("${prefix}_core_genes.tsv"), emit: gene_fam
    tuple val(ref_meta), path("${prefix}_feat_seqs/*.fasta"), emit: feat_seq, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${ref_meta.id}"
    """
    subset_core_gene.R $gene_fam $feat_seqs $sample_data $min_core_genes $min_core_samps $min_core_refs $max_core_genes ${prefix}_core_genes.tsv ${prefix}_feat_seqs
    """
}
