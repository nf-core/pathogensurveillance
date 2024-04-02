process R2TDIR {
    tag "$ref_meta.id"
    label 'process_single'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(ref_meta), path(busco_dir), path(ntFiles), path(aaFiles)

    output:
    tuple val(ref_meta), path("${prefix}/buscos"), emit: markers
    tuple val(ref_meta), path("${prefix}/dna_ref.fa"), emit: dna_ref
    
    //nfFiles = Channel.fromPath( '${busco_dir}/*ntformatted.fa')
    //aaFiles = Channel.fromPath( '${busco_dir}/*aaformatted.fa')

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${ref_meta.id}"
    """
    mkdir -p ${prefix}/buscos
    cat ${ntFiles.join(' ')} > ${prefix}/dna_ref.fa
    cp ${aaFiles.join(' '} > ${prefix}/buscos/
    """
}
