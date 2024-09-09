process BBMAP_SENDSKETCH {
    tag "$meta.id"
    label 'process_single'
    maxForks 1

    conda "bioconda::bbmap=39.08"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap:39.01--h5c4e2a8_0':
        'biocontainers/bbmap:39.01--h5c4e2a8_0' }"

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path('*.txt'), emit: hits
    tuple val(meta), env(DEPTH)   , emit: depth
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def file_used = file.size() > 1 ? file[0] : file
    """
    sendsketch.sh \\
        in=${file_used} \\
        out=${prefix}.txt \\
        $args \\
        -Xmx${task.memory.toGiga()}g

    DEPTH=\$(sed -n -e 's/.*Depth: \\([0-9.]\\+\\).*/\\1/p' ${prefix}.txt)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
}
