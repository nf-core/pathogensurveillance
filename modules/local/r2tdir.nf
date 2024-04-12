process R2TDIR {
    tag "all"
    label 'process_single'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    path busco_dir

    output:
    path "r2t_markers", emit: markers
    path "dna_ref.fa", emit: dna_ref 

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir r2t_markers
    cat **/*ntformatted.fa > dna_ref.fa
    cp **/*aaformatted.fa r2t_markers/
    """
}
