process TRIM_AND_SKETCH {
    tag "${meta.id}"
    label 'process_low'
    maxRetries 5
    maxErrors 30

    conda "bioconda::khmer=3.0.0a3 bioconda::sourmash=4.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/khmer:3.0.0a3--py37haa7609a_2' :
        'docker.io/zacharyfoster/khmer_sourmash:0.1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.sig"), emit: signatures
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    trim-low-abund.py \\
        -M ${task.memory.toGiga()}G \\
        $args \\
        $fasta

    sourmash sketch \\
        $args2 \\
        --merge '${prefix}' \\
        --output '${prefix}.sig' \\
        *.abundtrim

    rm *.abundtrim

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        khmer: \$( trim-low-abund.py --version 2>&1 | grep ^khmer | sed 's/^khmer //' )
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' )
    END_VERSIONS
    """
}

