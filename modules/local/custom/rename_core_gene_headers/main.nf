process RENAME_CORE_GENE_HEADERS {
    tag "$ref_meta.id"
    label 'process_single'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(ref_meta), path(feat_seqs)

    output:
    tuple val(ref_meta), path("${prefix}_feat_seqs_renamed"), emit: feat_seqs
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${ref_meta.id}"
    """
    # Create folder for output
    mkdir ${prefix}_feat_seqs_renamed

    # Rename headers to just sample ID
    for file in ${feat_seqs}/*.fasta
    do
        sed 's/>.*genome:\\(.*\\)gene:.*/>\\1/g' \$file > ${prefix}_feat_seqs_renamed/\$(basename \$file)
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(echo \$(sed --version | head -n 1 | awk '{print \$4}'))
    END_VERSIONS
    """
}
