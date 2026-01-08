process MULTIQC {
    tag "$meta.id"
    label 'process_single'

    // NOTE: hold at 1.28 to avoid "Illegal instruction" error on OSU CQLS HPC Slurm cluster
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c6c120d559d7ee04c7442b61ad7cf5a9e8970be5feefb37d68eeaa60c1034eb/data' :
        'community.wave.seqera.io/library/multiqc:1.32--d58f60e4deb769bf' }"

    input:
    tuple val(meta), path(multiqc_files, stageAs: "?/*")
    path(multiqc_config)
    path(extra_multiqc_config)
    path(multiqc_logo)
    path(replace_names)
    path(sample_names)

    output:
    tuple val(meta), path("${prefix}_multiqc")                     , emit: outdir
    tuple val(meta), path("${prefix}_multiqc/*multiqc_report.html"), emit: report
    tuple val(meta), path("${prefix}_multiqc/*_data")              , emit: data
    tuple val(meta), path("${prefix}_multiqc/*_plots")             , optional:true, emit: plots
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def config = multiqc_config ? "--config $multiqc_config" : ''
    def extra_config = extra_multiqc_config ? "--config $extra_multiqc_config" : ''
    def logo = multiqc_logo ? "--cl-config 'custom_logo: \"${multiqc_logo}\"'" : ''
    def replace = replace_names ? "--replace-names ${replace_names}" : ''
    def samples = sample_names ? "--sample-names ${sample_names}" : ''
    """
    multiqc \\
        --force \\
        --outdir ${prefix}_multiqc \\
        $args \\
        $config \\
        $extra_config \\
        $logo \\
        $replace \\
        $samples \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """

    stub:
    """
    mkdir multiqc_data
    mkdir multiqc_plots
    touch multiqc_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
