process VCF_TO_SNPALN {
    tag "$ref_meta.id"
    label 'process_low'

    conda "bioconda::perl-bioperl=1.7.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl-bioperl:1.7.8--hdfd78af_1':
        'biocontainers/perl-bioperl:1.7.8--hdfd78af_1' }"

    input:
    tuple val(ref_meta), path(tab)

    output:
    tuple val(ref_meta), path("*.fasta"), emit: fasta
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${ref_meta.id}"
    """
    vcftab_to_snpaln_nodel.pl --output_ref -i ${tab} > ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl --version | grep 'This is perl' | sed 's/.*(v//g' | sed 's/)//g')
    END_VERSIONS
    """
}
