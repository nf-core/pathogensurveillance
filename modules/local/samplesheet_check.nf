/*
Validates the input data and returns a reformatted version that is used for the rest of the pipeline.
*/

process SAMPLESHEET_CHECK {
    tag "$samplesheet"

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    path samplesheet

    output:
    path '*.csv'       , emit: csv
    path "versions.yml", emit: versions

    script:
    """
    check_samplesheet.py $samplesheet samplesheet.valid.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
