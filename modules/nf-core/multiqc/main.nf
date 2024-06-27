process MULTIQC {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.18--pyhdfd78af_0' :
        'biocontainers/multiqc:1.18--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(multiqc_files, stageAs: "?/*")
    path(multiqc_config)
    path(extra_multiqc_config)
    path(multiqc_logo)

    output:
    tuple val(meta), path("${prefix}_multiqc"), emit: outdir
    //tuple val(meta), path("*multiqc_report.html"), emit: report
    //tuple val(meta), path("*_data")              , emit: data
    //tuple val(meta), path("*_plots")             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def config = multiqc_config ? "--config $multiqc_config" : ''
    def extra_config = extra_multiqc_config ? "--config $extra_multiqc_config" : ''
    def logo = multiqc_logo ? /--cl-config 'custom_logo: "${multiqc_logo}"'/ : ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    multiqc \\
        --force \\
        --outdir ${prefix}_multiqc \\
        $args \\
        $config \\
        $extra_config \\
        $logo \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """

    stub:
    """
    touch multiqc_data
    touch multiqc_plots
    touch multiqc_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
