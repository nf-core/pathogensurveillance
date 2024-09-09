process VCF_TO_TAB {
    tag "$ref_meta.id"
    label 'process_low'

    conda "bioconda::vcftools=0.1.16"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcftools:0.1.16--pl5321h9a82719_6':
        'biocontainers/vcftools:0.1.16--pl5321h9a82719_6' }"

    input:
    tuple val(ref_meta), path(vcf)

    output:
    tuple val(ref_meta), path("*.tab"), emit: tab
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${ref_meta.id}"
    """
    gunzip -c ${vcf} | vcf-to-tab > ${prefix}.tab
    sed -i 's_/\\S\\+_/_g'  ${prefix}.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcftools: \$(echo \$(vcftools --version 2>&1) | sed 's/^.*VCFtools (//;s/).*//')
    END_VERSIONS
    """
}
