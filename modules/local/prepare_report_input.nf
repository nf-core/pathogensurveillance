process PREPARE_REPORT_INPUT {
    tag "$group_meta.id"
    label 'process_low'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(group_meta), val(ref_metas), file(sendsketchs), file(ref_data), file(quast_dirs), file(vcfs), file(snp_aligns), file(snp_phylos), file(ani_matrix), file(core_phylo), file(pocp), file(assigned_refs)
    path samp_data
    path multiqc_data
    path multiqc_plots
    path multiqc_report
    path versions
    path messages

    output:
    tuple val(group_meta), path("${prefix}_inputs"), emit: report_input

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${group_meta.id}"
    """
    # Make directory for ${prefix}_inputs so that a single path can be passed as parameters
    mkdir ${prefix}_inputs

    # Put multiqc's output into a single folder for organization
    mkdir ${prefix}_inputs/multiqc
    cp -r ${multiqc_data} ${prefix}_inputs/multiqc/
    cp -r ${multiqc_plots} ${prefix}_inputs/multiqc/
    cp -r ${multiqc_report} ${prefix}_inputs/multiqc/

    # Put quast's output into a single folder for organization
    mkdir ${prefix}_inputs/quast
    if [ ! -z "${quast_dirs}" ]; then
      cp -r ${quast_dirs} ${prefix}_inputs/quast/
    fi

    # Put sendsketch's output into a single folder for organization
    mkdir ${prefix}_inputs/sendsketch
    cp -r ${sendsketchs} ${prefix}_inputs/sendsketch/

    # Put variant data into a single folder for organization
    mkdir ${prefix}_inputs/variant_data
    if [ ! -z "${snp_phylos}" ]; then
         cp -r ${snp_phylos} ${prefix}_inputs/variant_data/
    fi
    if [ ! -z "${vcfs}" ]; then
        cp -r ${vcfs} ${prefix}_inputs/variant_data/
    fi
    if [ ! -z "${snp_aligns}" ]; then
        cp -r ${snp_aligns} ${prefix}_inputs/variant_data/
    fi

    # Save report group name to file
    echo "${group_meta.id}" > ${prefix}_inputs/group_id.txt

    # Move RefSeq reference data for each sample (for phylogenetic context) to their own directory
    mkdir ${prefix}_inputs/ref_data
    cp -r ${ref_data} ${prefix}_inputs/ref_data/

    # Move other single-value paths to input directory
    mkdir other_${prefix}_inputs
    mv ${samp_data} ${prefix}_inputs/samp_data.csv
    mv ${ani_matrix} ${prefix}_inputs/ani_matrix.csv
    mv ${assigned_refs} ${prefix}_inputs/assigned_refs.csv
    if [ ! -z "${core_phylo}" ]; then
        mv ${core_phylo} ${prefix}_inputs/core_phylo.treefile
    fi
    mv ${pocp} ${prefix}_inputs/pocp.tsv
    mv ${versions} ${prefix}_inputs/versions.yml
    mv ${messages} ${prefix}_inputs/messages.tsv
    """
}
