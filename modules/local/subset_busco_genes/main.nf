process SUBSET_BUSCO_GENES {
    tag "$group_meta.id"
    label 'process_single'

    conda "conda-forge::r-base=4.2.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'biocontainers/r-base:4.2.1' }"

    input:
    tuple val(group_meta), path(busco_out_dirs), val(samp_ref_pairs_raw)
    val min_genes
    val max_genes

    output:
    tuple val(group_meta), path("${prefix}_feat_seqs/${prefix}--cluster_*"), optional: true, emit: feat_seqs
    tuple val(group_meta), path("message_data.tsv")             , emit: message_data
    path "versions.yml"                                         , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${group_meta.id}"
    """
    echo "${samp_ref_pairs_raw}" | base64 -d > ${prefix}_samp_ref_pairs.tsv

    subset_busco_gene.R ${prefix}_samp_ref_pairs.tsv $min_genes $max_genes ${prefix} ${busco_out_dirs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
