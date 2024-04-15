process R2TBIN {
    tag "markers"
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::biopython=1.78" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'biocontainers/biopython:1.78' }"

    input:
    path markers

    output:
    path "busco_markers", emit: busco_markers 

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir busco_markers
    r2tbinder.py
    """
}
