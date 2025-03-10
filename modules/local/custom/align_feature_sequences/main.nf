process ALIGN_FEATURE_SEQUENCES {
    tag "$ref_meta.id"
    label 'process_low'

    conda "conda-forge::mafft=7.526 bioconda::perl-bioperl=1.7.8 conda-forge::parallel=20230522"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'zachary-foster/mafft-perl':
        'zachary-foster/mafft-perl' }"

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
        mafft: \$(mafft --version 2>&1 | sed 's/^v//' | sed 's/ (.*)//')
    END_VERSIONS
    """
}

