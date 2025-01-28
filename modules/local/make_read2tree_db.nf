process MAKE_READ2TREE_DB {
    tag "$report_meta.id"
    label 'process_low'

    conda "conda-forge::r-base=4.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'biocontainers/r-base:4.2.1' }"

    input:
    tuple val(report_meta), path(busco_dir)
    val lineage_db

    output:
    tuple val(report_meta), path("${prefix}_read2tree_db_dna.fa")   , emit: ref_dna
    tuple val(report_meta), path("${prefix}_read2tree_db_aa")       , emit: ref_aa
    tuple val(report_meta), path("${prefix}_read2tree_ref_meta.csv"), emit: ref_meta
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${report_meta.id}"
    """
    make_read2tree_db.R ${lineage_db} ${prefix}_read2tree_db_dna.fa ${prefix}_read2tree_db_aa ${prefix}_read2tree_ref_meta.csv ${busco_dir}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
