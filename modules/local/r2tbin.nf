process R2TBIN {
    tag "$markers"
    label 'process_single'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    path markers

    output:
    path "busco_markers", emit: busco_markers 

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir busco_markers
    python r2tbinder.py
    mv *.fasta busco_markers/
    """
}
