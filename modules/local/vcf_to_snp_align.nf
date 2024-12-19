process VCF_TO_SNP_ALIGN {
    tag "$ref_meta.id"
    label 'process_single'

    conda "conda-forge::r-base=4.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'biocontainers/r-base:4.2.1' }"

    input:
    tuple val(ref_meta), path(vcf)

    output:
    tuple val(ref_meta), path("${prefix}.fasta")       , emit: fasta
    tuple val(ref_meta), path("removed_sample_ids.txt"), emit: removed_sample_ids
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${ref_meta.id}"
    """
    vcf_to_snp_align.R ${vcf} ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
