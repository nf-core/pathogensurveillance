process PARSE_ASSEMBLIES {
    tag "$taxon"
    label 'process_single'

    conda "conda-forge::r-rcppsimdjson=0.1.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3c/3c9b7a9283feb72d16c22aa52a48404bc346769643c16d3f09830cd4955e89cf/data':
        'community.wave.seqera.io/library/r-rcppsimdjson:0.1.12--e12a2b75de86869a' }"

    container "quay.io/nf-core/rcppsimdjson:0.2"

    input:
    tuple val(taxon), path(json)

    output:
    tuple val(taxon), path("${prefix}.tsv"), emit: stats
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: taxon
    """
    parse_assemblies.R ${json} ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        rcppsimdjson: \$(echo \$(Rscript -e "cat(format(packageVersion('RcppSimdJson')))"))
    END_VERSIONS
    """
}
