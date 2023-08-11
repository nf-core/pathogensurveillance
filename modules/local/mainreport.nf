process MAIN_REPORT {
    tag "$group_meta.id"
    label 'process_low'

    conda null // TODO: it just uses the local computers R for now
    container null

    input:
    tuple val(group_meta), val(ref_metas), file(snp_phylos), file(ani_matrix), file(core_phylo)
    path samp_data
    path ref_data

    output:
    tuple val(group_meta), path("main_report/${prefix}_report"), emit: html
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${group_meta.id}"
    """
    # Copy source of report here cause bookdown seems to want to make its output in the source
    cp -r ${projectDir}/assets/main_report ./

    Rscript -e "workdir<-getwd()
        bookdown::render_book('main_report',
        output_format = 'all',
        params = list(
        samp_data = \\\"$samp_data\\\",
        ref_data = \\\"$ref_data\\\",
        snp_phylos = \\\"$snp_phylos\\\",
        ani_matrix = \\\"$ani_matrix\\\",
        core_phylo = \\\"$core_phylo\\\"
        ),
        output_dir = \\\"${prefix}_report\\\",
        knit_root_dir = workdir)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rmarkdown: \$(Rscript -e "cat(paste(packageVersion('rmarkdown'), collapse='.'))")
    END_VERSIONS
    """
}
