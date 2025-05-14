process PICARD_FORMAT {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::picard=3.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.1.1--hdfd78af_0' :
        'biocontainers/picard:3.1.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(reference), path(ref_index)

    output:
    tuple val(meta), path("${prefix}.bam")        , emit: bam
    tuple val(meta), path("${prefix}.bai")        , emit: bai,     optional: true
    tuple val(meta), path("${prefix}.MarkDuplicates.metrics.txt"), emit: metrics
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    prefix = task.ext.prefix    ?: "${meta.id}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    """
    picard \\
        -Xmx${avail_mem}M \\
        AddOrReplaceReadGroups \\
        $args \\
        --INPUT ${bam} \\
        --OUTPUT ${prefix}_AddOrReplaceReadGroups.bam

    picard \\
        SortSam \\
        -Xmx${avail_mem}M \\
        --INPUT ${prefix}_AddOrReplaceReadGroups.bam \\
        --OUTPUT ${prefix}_SortSam_1.bam \\
        --SORT_ORDER coordinate \\
        $args2

    picard \\
        -Xmx${avail_mem}M \\
        MarkDuplicates \\
        $args3 \\
        --INPUT ${prefix}_SortSam_1.bam \\
        --OUTPUT ${prefix}_MarkDuplicates.bam \\
        --REFERENCE_SEQUENCE $reference \\
        --METRICS_FILE ${prefix}.MarkDuplicates.metrics.txt

    picard \\
        SortSam \\
        -Xmx${avail_mem}M \\
        --INPUT ${prefix}_MarkDuplicates.bam \\
        --OUTPUT ${prefix}.bam \\
        --SORT_ORDER coordinate \\
        $args2

    rm ${prefix}_AddOrReplaceReadGroups.bam
    rm ${prefix}_SortSam_1.bam
    rm ${prefix}_MarkDuplicates.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard AddOrReplaceReadGroups --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix    ?: "${meta.id}"

    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard AddOrReplaceReadGroups --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
