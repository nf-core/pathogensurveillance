/*
Validates the input data and returns a reformatted version that is used for the rest of the pipeline.
*/

process SAMPLESHEET_CHECK {
    tag "$samplesheet"

    conda "conda-forge::r-rentrez=1.2.3"                                           
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :            
        'quay.io/biocontainers/r-entrez:1.1.0--r3.4.1_0' }"                                 

    input:
    path sample_csv
    path reference_csv

    output:
    path 'sample_metadata.csv'    , emit: sample_metadata
    path 'reference_metadata.csv' , emit: reference_metadata
    path "versions.yml"           , emit: versions

    script:
    """
    check_samplesheet.R ${sample_metadata} ${reference_metadata} sample_metadata.csv reference_metadata.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
