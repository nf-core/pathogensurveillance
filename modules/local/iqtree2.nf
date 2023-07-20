process IQTREE2 {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::iqtree=2.1.4_beta"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/iqtree:2.1.4_beta--hdcc8f71_0' :
        'biocontainers/iqtree:2.1.4_beta--hdcc8f71_0' }"

    input:
    tuple val(meta), path(alignment, stageAs: 'input_alignments/*')
    val constant_sites

    output:
    tuple val(meta), path("*.treefile"), emit: phylogeny
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def fconst_args = constant_sites ? "-fconst $constant_sites" : ''
    def memory      = task.memory.toString().replaceAll(' ', '').replaceAll('B', '')
    def first_align = alignment[0]
    """
    # Get number of samples to decide whether or not to bootstrap
    NSAMPLE=\$(grep '>' $first_align | wc -l)
    if [[\$NSAMPLE -gt 3]]; then
        BOOT="-B 1000"
    else
        BOOT=""
    fi

    # Create phylogenetic tree
    iqtree2 \\
        $fconst_args \\
        \$BOOT \\
        $args \\
        -s input_alignments \\
        -nt AUTO \\
        -ntmax $task.cpus \\
        -mem $memory \\

    # Save version information
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(echo \$(iqtree -version 2>&1) | sed 's/^IQ-TREE multicore version //;s/ .*//')
    END_VERSIONS
    """
}
