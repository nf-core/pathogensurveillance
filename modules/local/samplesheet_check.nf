/*
Validates the input data and returns a reformatted version that is used for the rest of the pipeline.
*/

process SAMPLESHEET_CHECK {
    tag "$samplesheet"

    conda "conda-forge::r-base=4.2.1"                                           
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :            
        'quay.io/biocontainers/r-base:4.2.1' }"                                 

    input:
    path samplesheet

    output:
    path '*.csv'       , emit: csv
    path "versions.yml", emit: versions

    script:
    """
    check_samplesheet.R $samplesheet samplesheet.valid.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
