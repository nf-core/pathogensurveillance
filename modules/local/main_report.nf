process MAIN_REPORT {
    tag "$group_meta.id"
    label 'process_low'

    conda "bioconda::bioconductor-ggtree=3.10.0 conda-forge::quarto=1.5.55 conda-forge::r-ade4=1.7_22 conda-forge::r-adegenet=2.1.10 conda-forge::r-ape=5.8 conda-forge::r-biocmanager=1.30.25 conda-forge::r-devtools=2.4.5 conda-forge::r-dplyr=1.1.4 conda-forge::r-ggdendro=0.2.0 conda-forge::r-ggnewscale=0.5.0 conda-forge::r-ggplot2=3.5.1 conda-forge::r-heatmaply=1.5.0 conda-forge::r-igraph=2.0.3 conda-forge::r-kableextra=1.4.0 conda-forge::r-knitr=1.48 conda-forge::r-leaflet=2.2.2 conda-forge::r-magrittr=2.0.3 conda-forge::r-palmerpenguins=0.1.1 conda-forge::r-phangorn=2.11.1 conda-forge::r-pheatmap=1.0.12 conda-forge::r-plotly=4.10.4 conda-forge::r-purrr=1.0.2 conda-forge::r-rcrossref=1.2.0 conda-forge::r-readr=2.1.5 conda-forge::r-tidyverse=2.0.0 conda-forge::r-visnetwork=2.1.2 conda-forge::r-webshot2=0.1.1 conda-forge::r-yaml=2.3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/main-report-r-packages':
        'docker.io/zacharyfoster/main-report-r-packages:0.13' }"

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
    # Copy source of report here cause quarto seems to want to make its output in the source
    cp -r --dereference main_report_template main_report

    # Render the report
    quarto render main_report \\
        ${args} \\
        --output-dir ${prefix}_report \\
        -P inputs:../${inputs}

    # Rename outputs
    mv main_report/${prefix}_report/index.html ${prefix}_pathsurveil_report.html
    #mv main_report/${prefix}_report/index.pdf ${prefix}_pathsurveil_report.pdf

    # Save version of quarto used
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto --version)
    END_VERSIONS
    """
}
