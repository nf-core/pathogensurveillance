process PREPARE_REPORT_INPUT {
    tag "$group_meta.id"
    label 'process_low'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(group_meta), path(sample_data), path(ref_data), path(sendsketch), path(ncbi_ref_meta), path(selected_refs), path(ani_matrix), path(mapping_ref), path(snp_aligns), path(snp_phylos), path(core_phylo_refs, stageAs: 'core_phylo_refs.csv'), path(pocp), path(core_phylos, stageAs: 'core_phylos/*'), path(busco_refs, stageAs: 'busco_refs.csv'), path(busco_phylos, stageAs: 'busco_phylos/*'), path(multiqc), path(messages)
    path versions

    output:
    tuple val(group_meta), path("${prefix}_inputs"), emit: report_input

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${group_meta.id}"
    config_file_text = workflow.configFiles.collect{"  - ${it}"}.join("\n")
    """
    # Make directory for ${prefix}_inputs so that a single path can be passed as parameters
    mkdir ${prefix}_inputs

    # Add sample metadata for this report group
    cp -r ${sample_data} ${prefix}_inputs/sample_data.csv

    # Add reference metadata for this report group
    cp -r ${ref_data} ${prefix}_inputs/reference_data.csv

    # Put sendsketch's output into a single folder for organization
    mkdir ${prefix}_inputs/sendsketch
    cp -r ${sendsketch} ${prefix}_inputs/sendsketch/

    # Put all of the statistics for NCBI references considered into a single folder
    if [ ! -z "${ncbi_ref_meta}" ]; then
        mkdir ${prefix}_inputs/ncbi_reference_data
        cp -r ${ncbi_ref_meta} ${prefix}_inputs/ncbi_reference_data/
    fi

    # Put the metadata for references selected for each sample into a single folder
    if [ ! -z "${selected_refs}" ]; then
        mkdir ${prefix}_inputs/selected_references
        cp -r ${selected_refs} ${prefix}_inputs/selected_references/
    fi

    # Add estimated ANI matrix from sourmash
    mv ${ani_matrix} ${prefix}_inputs/sourmash_ani_matrix.csv

    # Add metadata for references assined for variant calling
    if [ ! -z "${mapping_ref}" ]; then
        mv ${mapping_ref} ${prefix}_inputs/mapping_references.csv
    fi

    # Add SNP alignment from variant calling
    if [ ! -z "${snp_aligns}" ]; then
        mkdir ${prefix}_inputs/snp_alignments
        cp -r ${snp_aligns} ${prefix}_inputs/snp_alignments/
    fi

    # Add SNP phylogenies from variant calling
    if [ ! -z "${snp_phylos}" ]; then
        mkdir ${prefix}_inputs/snp_trees
        cp -r ${snp_phylos} ${prefix}_inputs/snp_trees/
    fi

    # Add seleted references for the core gene phylogeny
    if [ ! -z "${core_phylo_refs}" ]; then
        cp -r ${core_phylo_refs} ${prefix}_inputs/core_gene_tree_references.csv
    fi

    # Add POCP estimate from the core genome phylogeny
    if [ ! -z "${pocp}" ]; then
        cp -r ${pocp} ${prefix}_inputs/pocp.csv
    fi

    # Add core genome phylogenies
    if [ ! -z "${core_phylos}" ]; then
        mkdir ${prefix}_inputs/core_gene_trees
        cp -r ${core_phylos} ${prefix}_inputs/core_gene_trees/
    fi

    # Add seleted references for the busco phylogeny
    if [ ! -z "${busco_refs}" ]; then
        cp -r ${busco_refs} ${prefix}_inputs/busco_tree_references.csv
    fi

    # Add busco phylogeny
    if [ ! -z "${busco_phylos}" ]; then
        cp -r ${busco_phylos} ${prefix}_inputs/busco_tree.nwk
    fi

    # Put multiqc's output into a single folder for organization
    cp -r ${multiqc} ${prefix}_inputs/multiqc

    # Add pipeline status messages
    if [ ! -z "${messages}" ]; then
        mv ${messages} ${prefix}_inputs/messages.csv
    fi

    # Add versions of software used
    mv ${versions} ${prefix}_inputs/versions.yml

    # Record run metadata
    cat <<-END_INFO > ${prefix}_inputs/pathogensurveillance_run_info.yml
    group_id: ${group_meta.id}
    command_line: ${workflow.commandLine}
    commit_id: ${workflow.commitId}
    config_files:
    ${config_file_text}
    container_engine: ${workflow.containerEngine}
    profile: ${workflow.profile}
    revision: ${workflow.revision}
    run_name: ${workflow.runName}
    session_id: ${workflow.sessionId}
    start_time: ${workflow.start}
    nextflow_version: ${nextflow.version}
    pipeine_verison: ${workflow.manifest.version}
    END_INFO
    """
}
