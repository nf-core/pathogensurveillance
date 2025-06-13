process MAIN_REPORT {
    tag "$group_meta.id"
    label 'process_low'

    conda "conda-forge::quarto=1.6.41 bioconda::r-pathosurveilr=0.4.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/zacharyfoster/main-report-r-packages:0.20':
        'docker.io/zacharyfoster/main-report-r-packages:0.20' }"

    input:
    tuple val(group_meta), file(inputs)
    path template, stageAs: 'main_report_template'

    output:
    tuple val(group_meta), path("${prefix}_pathsurveil_report.html"), emit: html
    tuple val(group_meta), path("${prefix}_pathsurveil_report.pdf") , emit: pdf, optional: true
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${group_meta.id}"
    """
    ls -la ${inputs}/
    ls -la ${inputs}/trees/

    # Tell quarto where to put cache so it does not try to put it where it does not have permissions
    export XDG_CACHE_HOME="\$(pwd)/cache"

    # Copy source of report here cause quarto seems to want to make its output in the source
    cp -r --dereference main_report_template main_report

    # Render the report
    quarto render main_report \\
        ${args} \\
        --output-dir ${prefix}_report \\
        -P inputs:../${inputs}

    # Rename outputs
    mv main_report/${prefix}_report/index.html ${prefix}_pathsurveil_report.html

    # Clean up
    rm -r main_report

    # Save version of quarto used
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto --version)
        r-PathoSurveilR: \$(Rscript -e "cat(as.character(packageVersion('dplyr')))")
    END_VERSIONS
    """
}
