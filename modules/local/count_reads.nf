process COUNT_READS {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(fastqs)

    output:
    tuple val(meta), env(READ_COUNT), emit: read_count
    tuple val(meta), env(READ_LEN), emit: total_bp
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    READ_COUNT=\$(zgrep -c '@' ${fastqs[0]})
    READ_LEN=\$(zgrep -m 1 '^[^@+#]' ${fastqs[0]} | wc -m)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        zgrep: \$(zgrep --version 2>&1) | sed 's/^.*gzip) //; s/ .*\$//'
    END_VERSIONS
    """
}

