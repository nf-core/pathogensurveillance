/*
Validates the input data and returns a reformatted version that is used for the rest of the pipeline.
*/

process SAMPLESHEET_CHECK {
    tag "$sample_tsv"

    conda "conda-forge::quarto=1.6.41 bioconda::r-pathosurveilr=0.4.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/zacharyfoster/main-report-r-packages:0.24':
        'docker.io/zacharyfoster/main-report-r-packages:0.24' }"

    input:
    path sample_tsv
    path reference_tsv
    val max_samples

    output:
    path 'sample_metadata.tsv'   , emit: sample_data
    path 'reference_metadata.tsv', emit: reference_data
    path 'message_data.tsv'      , emit: message_data
    path "versions.yml"          , emit: versions

    script:
    def entrez_key_set = secrets.NCBI_API_KEY ? "export ENTREZ_KEY='${secrets.NCBI_API_KEY}'" : ''
    """
    ${entrez_key_set}

    check_samplesheet.R ${sample_tsv} ${max_samples} ${reference_tsv}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-PathoSurveilR: \$(Rscript -e "cat(as.character(packageVersion('dplyr')))")
    END_VERSIONS
    """
}
