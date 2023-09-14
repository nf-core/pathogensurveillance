process MAIN_REPORT {
    tag "$group_meta.id"
    label 'process_low'

    conda "conda-forge::quarto=1.3.450" // TODO: it just uses the local computers R packages for now
    container null

    input:
    tuple val(group_meta), val(ref_metas), file(snp_phylos), file(snp_align), file(ani_matrix), file(core_phylo)
    path samp_data
    path ref_data

    output:
    tuple val(group_meta), path("main_report/${prefix}_report"), emit: html
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${group_meta.id}"
    """
    # Copy source of report here cause quarto seems to want to make its output in the source
    cp -r ${projectDir}/assets/main_report ./
    
    quarto render main_report \\
        --output-dir ${prefix}_report \\
        -P samp_data:../${samp_data} \\
        -P ref_data:../${ref_data} \\
        -P snp_phylos:../${snp_phylos} \\
        -P ani_matrix:../${ani_matrix} \\
        -P core_phylo:../${core_phylo} \\
        -P snp_align:../${snp_align}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto --version)
    END_VERSIONS
    """
}
