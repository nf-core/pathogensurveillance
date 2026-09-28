/*
Validates the input data and returns a reformatted version that is used for the rest of the pipeline.
*/

process SAMPLESHEET_CHECK {
    tag "input metadata"

    conda "conda-forge::quarto=1.6.41 bioconda::r-pathosurveilr=0.4.8"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e4/e487168aaa5f7b7a2dfabc1f308869fac1c466ab6bd34acd7a3715d927337c74/data':
        'community.wave.seqera.io/library/r-pathosurveilr_quarto:d4f39be8e8ae4734' }"

    input:
    path sample_tsv    , stageAs: 'input_sample_metadata.txt'
    path reference_tsv , stageAs: 'input_reference_metadata.txt'
    path report_tsv    , stageAs: 'input_report_metadata.txt'
    val max_samples

    output:
    path 'sample_metadata.tsv'   , emit: sample_data
    path 'reference_metadata.tsv', emit: reference_data
    path 'report_metadata.tsv'   , emit: report_data
    path 'message_data.tsv'      , emit: message_data
    path "versions.yml"          , emit: versions_samplesheet_check, topic: versions

    script:
    def entrez_key_set = secrets.NCBI_API_KEY ? "export ENTREZ_KEY='${secrets.NCBI_API_KEY}'" : ''
    """
    ${entrez_key_set}

    check_samplesheet.R input_sample_metadata.txt ${max_samples} input_reference_metadata.txt input_report_metadata.txt ${projectDir}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-PathoSurveilR: \$(Rscript -e "cat(as.character(packageVersion('dplyr')))")
    END_VERSIONS
    """
}
