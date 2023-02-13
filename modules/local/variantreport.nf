process VARIANTREPORT {
    tag "$ref_meta.id"
    label 'process_low'

    conda "conda-forge::r-base=4.1.0 conda-forge::r-rmarkdown=2.9 conda-forge::r-yaml=2.2.1 bioconda::r-poppr=2.8.5 bioconda::r-adegenet=2.1.1 r::r-ape=5.6_2 r::r-tidyverse=1.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-31ad840d814d356e5f98030a4ee308a16db64ec5:0e852a1e4063fdcbe3f254ac2c7469747a60e361-0' :
        'quay.io/biocontainers/mulled-v2-31ad840d814d356e5f98030a4ee308a16db64ec5:0e852a1e4063fdcbe3f254ac2c7469747a60e361-0' }"

    input:
    tuple val(ref_meta), path(fasta)
    path samplesheet

    output:
    tuple val(ref_meta), path("*.html"), emit: html
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${ref_meta.id}"
    """
    Rscript -e "workdir<-getwd()
        rmarkdown::render('$projectDir/assets/test_report.Rmd',
        params = list(
        fasta_align = \\\"$fasta\\\",
        metadata = \\\"$samplesheet\\\"
        ),
        knit_root_dir=workdir,
        output_dir=workdir)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rmarkdown: \$(Rscript -e "cat(paste(packageVersion('rmarkdown'), collapse='.'))")
    END_VERSIONS
    """
}
