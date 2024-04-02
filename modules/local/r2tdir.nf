process R2TDIR {
    tag "$ref_meta.id"
    label 'process_single'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(ref_meta), path(busco_dir)

    output:
    tuple val(ref_meta), path("buscos"), emit: markers
    tuple val(ref_meta), path("dna_ref.fa"), emit: dna_ref 

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${ref_meta.id}"
    """
    mkdir -p ${busco_dir}/buscos
    cat ${busco_dir}/ > ${busco_dir}/dna_ref.fa
    cp ${busco_dir}/ > ${busco_dir}/buscos/
    """
}
