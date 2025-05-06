process EXTRACT_FEATURE_SEQUENCES {
    tag "$ref_meta.id"
    label 'process_low'

    conda "bioconda::pirate=1.0.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pirate:1.0.5--hdfd78af_0' :
        'biocontainers/pirate:1.0.5--hdfd78af_0' }"

    input:
    tuple val(ref_meta), path(pirate_results)

    output:
    tuple val(ref_meta), path("${prefix}_feature_sequences"), emit: feat_seqs
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${ref_meta.id}"
    """
    # Set dosage ultra high to include all high copy sequences per strain, just get everything
    align_feature_sequences_mod.pl --dosage 1000000 --full-annot -i PIRATE.gene_families.ordered.tsv -g modified_gffs/ -o ${prefix}_feature_sequences/ -p ${task.cpus} -n --align-off

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pirate: \$( echo \$( PIRATE --version 2>&1) | sed 's/PIRATE //' )
    END_VERSIONS
    """
}

