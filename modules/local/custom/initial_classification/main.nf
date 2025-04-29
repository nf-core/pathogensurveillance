process INITIAL_CLASSIFICATION {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::r-pathosurveilr=0.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-pathosurveilr:0.3.1--r44hdfd78af_0' :
        'quay.io/biocontainers/r-pathosurveilr:0.3.1--r44hdfd78af_0' }"

    input:
    tuple val(meta), path(hits)

    output:
    tuple val(meta), path("${prefix}_taxa_found.tsv"), emit: taxa_found
    tuple val(meta), path("${prefix}_taxon_data.tsv"), emit: taxon_data
    tuple val(meta), env(DOMAIN)                     , emit: domain
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    initial_classification.R $hits

    DOMAIN="\$(cat domain.txt)"

    mv taxa_found.tsv ${prefix}_taxa_found.tsv
    mv taxon_data.tsv ${prefix}_taxon_data.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        rentrez: \$(echo \$(Rscript -e "cat(format(packageVersion('rentrez')))"))
    END_VERSIONS
    """
}
