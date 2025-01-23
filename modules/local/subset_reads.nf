process SUBSET_READS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::seqkit=2.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.2.0--h9ee0642_0':
        'biocontainers/seqkit:2.2.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(fastqs), val(read_count)

    output:
    tuple val(meta), path("*_subset.fastq.gz"), emit: reads
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    for f in ${fastqs.join(' ')}
    do
        seqkit head ${args} --threads $task.cpus -n ${read_count} -o "\$(basename \$f .fastq.gz)_subset.fastq.gz" \$f
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
    END_VERSIONS
    """
}
