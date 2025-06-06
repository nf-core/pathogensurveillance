process PICK_ASSEMBLIES {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::r-pathosurveilr=0.4.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-pathosurveilr:0.3.1--r44hdfd78af_0' :
        'quay.io/biocontainers/r-pathosurveilr:0.3.1--r44hdfd78af_0' }"

    input:
    tuple val(meta), path(found_taxa), path(assem_data_tsvs)
    val n_ref_strains
    val n_ref_species
    val n_ref_genera
    val only_latin_binomial_refs

    output:
    tuple val(meta), path("${prefix}_formatted.tsv"), emit: formatted
    tuple val(meta), path("${prefix}.tsv")          , emit: metadata
    tuple val(meta), env(COUNT)                     , emit: line_count
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    pick_assemblies.R ${found_taxa} ${n_ref_strains} ${n_ref_species} ${n_ref_genera} ${only_latin_binomial_refs} ${prefix} ${assem_data_tsvs}
    COUNT=\$(cat ${prefix}.tsv | wc -l)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-PathoSurveilR: \$(Rscript -e "cat(as.character(packageVersion('dplyr')))")
    END_VERSIONS
    """
}
