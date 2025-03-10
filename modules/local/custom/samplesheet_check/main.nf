/*
Validates the input data and returns a reformatted version that is used for the rest of the pipeline.
*/

process SAMPLESHEET_CHECK {
    tag "$sample_tsv"

    conda "conda-forge::r-rentrez=1.2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/zacharyfoster/rentrez:0.1' :
        'docker.io/zacharyfoster/rentrez:0.1' }"

    input:
    path sample_tsv
    path reference_tsv

    output:
    path 'sample_metadata.tsv'   , emit: sample_data
    path 'reference_metadata.tsv', emit: reference_data
    path 'message_data.tsv'      , emit: message_data
    path "versions.yml"          , emit: versions

    script:
    """
    check_samplesheet.R ${sample_tsv} ${reference_tsv}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
