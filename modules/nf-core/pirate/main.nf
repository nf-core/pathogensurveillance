process PIRATE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::pirate=1.0.5 bioconda::perl-bioperl=1.7.8 bioconda::mcl=14.137"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pirate:1.0.4--hdfd78af_2' :
        'quay.io/biocontainers/pirate:1.0.5--hdfd78af_0' }"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("results/*")                                   , emit: results
    tuple val(meta), path("results/core_alignment.fasta"), optional: true, emit: aln
    path "versions.yml"                                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Rename .gff3 to .gff if needed
    find -regex .*\\.gff3\$ | sed 's/3\$//' | xargs -I {} mv {}3 {}

    # Run pirate on all .gff in input directory
    PIRATE \\
        $args \\
        --threads $task.cpus \\
        --input ./ \\
        --output results/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pirate: \$( echo \$( PIRATE --version 2>&1) | sed 's/PIRATE //' )
    END_VERSIONS
    """
}
