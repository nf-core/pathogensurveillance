process SUBSETCOREGENES {
    tag "$ref_meta.id"
    label 'process_single'
                                                                                
    conda "conda-forge::coreutils=9.1"                                          
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :            
        'ubuntu:20.04' }"                                                       

    input:
    tuple val(ref_meta), path(gene_fam), path(feat_seqs)

    output:
    tuple val(ref_meta), path("${prefix}_core_genes.tsv"), emit: gene_fam
    tuple val(ref_meta), path("${prefix}_feat_seqs/*.fasta"), emit: feat_seq, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${ref_meta.id}"
    """
    # Find number of samples based on the number of extra per-sample columns
    NCOL=\$(awk -F '\\t' '{print NF - 22; exit}' $gene_fam)
    
    # Remove rows representing gene families not present in all samples or that have paralogs
    awk -F '\\t' "\\\$7 == \$NCOL && \\\$10 == 1 && \\\$11 == 0 && \\\$12 == 0" $gene_fam > ${prefix}_core_genes.tsv

    # Make a directory with links to the sequences for genes that are single-copy and core
    mkdir ${prefix}_feat_seqs
    cut -f2 ${prefix}_core_genes.tsv | sed 's/\$/.fasta/' | xargs -I {} cp -s \$PWD/$feat_seqs/{} \$PWD/${prefix}_feat_seqs/{}
    """
}
