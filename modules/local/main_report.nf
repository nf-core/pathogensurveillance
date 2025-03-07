process MAIN_REPORT {
    tag "$group_meta.id"
    label 'process_low'

    conda "conda-forge::r-rcurl conda-forge::r-xml conda-forge::r-devtools conda-forge::r-leaflet conda-forge::r-biocmanager conda-forge::quarto conda-forge::r-knitr conda-forge::r-dplyr conda-forge::r-ggplot2 conda-forge::r-readr conda-forge::r-purrr conda-forge::r-yaml conda-forge::r-ape conda-forge::r-magrittr conda-forge::r-pheatmap conda-forge::r-heatmaply conda-forge::r-tidyverse conda-forge::r-palmerpenguins conda-forge::r-ade4 conda-forge::r-adegenet bioconda::bioconductor-ggtree conda-forge::r-igraph conda-forge::r-visnetwork conda-forge::r-phangorn conda-forge::r-ggnewscale conda-forge::r-kableextra conda-forge::r-plotly conda-forge::r-webshot2 conda-forge::r-ggdendro conda-forge::r-rcrossref"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/zacharyfoster/main-report-r-packages:0.18':
        'docker.io/zacharyfoster/main-report-r-packages:0.18' }"
        //singularity not currently supported. First run command [singularity build main-report-r-packages.sif docker://zacharyfoster/main-report-r-packages]

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
    #mv main_report/${prefix}_report/index.pdf ${prefix}_pathsurveil_report.pdf

    # Save version of quarto used
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto --version)
    END_VERSIONS
    """
}
