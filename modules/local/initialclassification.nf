process INITIALCLASSIFICATION {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(hits), path(reads), val(ref_meta), path(reference)

    output:
    tuple val(meta), env(TAXON), path(hits), path(reads), val(ref_meta), path(reference), emit: result
    tuple val(meta), env(CLASS),                                                          emit: classification
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    initial_classification.py $hits taxon.txt class.txt
    TAXON="\$(cat taxon.txt)"
    CLASS="\$(cat class.txt)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
