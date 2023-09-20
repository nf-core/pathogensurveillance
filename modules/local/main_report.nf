process MAIN_REPORT {
    tag "$group_meta.id"
    label 'process_low'

    conda "conda-forge::quarto=1.3.450" // TODO: it just uses the local computers R packages for now
    container null

    input:
    tuple val(group_meta), val(ref_metas), file(snp_phylos), file(snp_aligns), file(vcfs), file(quast_dirs), file(ani_matrix), file(core_phylo), file(sendsketchs)
    path samp_data
    path ref_data
    path multiqc_data
    path multiqc_plots
    path multiqc_report
    path versions

    output:
    tuple val(group_meta), path("main_report/${prefix}_report"), emit: html
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${group_meta.id}"
    def ref_ids = ref_metas.collect{ it.id }.join(';')
    """
    # Copy source of report here cause quarto seems to want to make its output in the source
    cp -r ${projectDir}/assets/main_report ./
    
    # Put multiqc's output into a single folder for organization
    mkdir multiqc
    cp -r ${multiqc_data} multiqc/
    cp -r ${multiqc_plots} multiqc/                                              
    cp -r ${multiqc_report} multiqc/
    
    # Put quast's output into a single folder for organization
    mkdir quast
    cp -r ${quast_dirs} quast/
    
    # Put sendsketch's output into a single folder for organization
    mkdir sendsketch
    cp -r ${sendsketchs} sendsketch/

    # Put variant data into a single folder for organization
    mkdir variant_data
    cp -r ${snp_phylos} variant_data/
    cp -r ${vcfs} variant_data/
    cp -r ${snp_aligns} variant_data/
    
    # Render the report
    quarto render main_report \\
        --output-dir ${prefix}_report \\
        -P group:${group_meta.id} \\
        -P refs:${ref_ids} \\
        -P samp_data:../${samp_data} \\
        -P ref_data:../${ref_data} \\
        -P sendsketch:../sendsketch
        -P variant_data:../variant_data \\
        -P ani_matrix:../${ani_matrix} \\
        -P core_phylo:../${core_phylo} \\
        -P multiqc:../multiqc \\
        -P quast:../quast \\
        -P versions:../versions

    # Save version of quarto used
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto --version)
    END_VERSIONS
    """
}
