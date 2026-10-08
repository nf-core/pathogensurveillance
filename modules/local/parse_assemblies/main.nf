process PARSE_ASSEMBLIES {
    tag "$taxon"
    label 'process_single'

    conda "conda-forge::r-rcppsimdjson=0.1.12"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1e/1e78da6657c8abaa07c228d7528b85e841d5f4bc9d0b276d7eef658a8a85a819/data' :
        'community.wave.seqera.io/library/r-rcppsimdjson:0.1.12--61f5cb2fd0b45fdd' }"

    input:
    tuple val(taxon), path(json)

    output:
    tuple val(taxon), path("${prefix}.tsv"), emit: stats
    path "versions.yml"                    , emit: versions_parse_assemblies, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: taxon
    """
    parse_assemblies.R ${json} ${prefix}.tsv
    ls -hl .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        rcppsimdjson: \$(echo \$(Rscript -e "cat(format(packageVersion('RcppSimdJson')))"))
    END_VERSIONS
    """
}
