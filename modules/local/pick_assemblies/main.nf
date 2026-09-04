process PICK_ASSEMBLIES {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::quarto=1.6.41 bioconda::r-pathosurveilr=0.4.7"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ec/ec0e2ca110b9875eab3997bac2682b5735bcbff78a6025b98816dd016685bdd4/data':
        'community.wave.seqera.io/library/r-pathosurveilr_quarto:dc06886ee7b6ddcd' }"

    input:
    tuple val(meta), path(found_taxa), path(child_taxa), path(assem_data_tsvs), val(excluded_accessions)
    val n_ref_strains
    val n_ref_species
    val n_ref_genera
    val only_latin_binomial_refs

    output:
    tuple val(meta), path("${prefix}_formatted.tsv"), emit: formatted
    tuple val(meta), path("${prefix}.tsv")          , emit: metadata
    tuple val(meta), env('COUNT')                   , emit: line_count
    path "versions.yml"                             , emit: versions_pick_assemblies, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def excluded_arg = excluded_accessions != 'NONE' ? "\"${excluded_accessions}\"" : 'NONE'
    """
    pick_assemblies.R ${found_taxa} ${child_taxa} ${n_ref_strains} ${n_ref_species} ${n_ref_genera} ${only_latin_binomial_refs} ${prefix} ${excluded_arg} ${assem_data_tsvs}
    COUNT=\$(cat ${prefix}.tsv | wc -l)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-PathoSurveilR: \$(Rscript -e "cat(as.character(packageVersion('dplyr')))")
    END_VERSIONS
    """
}
