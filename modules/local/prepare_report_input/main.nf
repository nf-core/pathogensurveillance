process PREPARE_REPORT_INPUT {
    tag "$group_meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data' :
        'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(group_meta), path(sample_data), path(ref_data), path(sendsketch), path(ncbi_ref_meta), path(selected_refs), path(ani_matrix), path(mapping_ref), path(snp_aligns), path(snp_phylos), path(core_phylo_refs, stageAs: 'core_phylo_refs.tsv'), path(pocp), path(core_phylos, stageAs: 'core_phylos/*'), path(busco_refs, stageAs: 'busco_refs.tsv'), path(busco_phylos, stageAs: 'busco_phylos/*'), path(multiqc), path(messages), path(versions)
    path output_format_json

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

    # Add output format metadata file for use with PathoSurveilR
    cp ${output_format_json} ${prefix}_inputs/

    # Add sample metadata for this report group
    mkdir -p ${prefix}_inputs/metadata
    cp -r ${sample_data} ${prefix}_inputs/metadata/sample_metadata.tsv

    # Add reference metadata for this report group
    mkdir -p ${prefix}_inputs/metadata
    cp -r ${ref_data} ${prefix}_inputs/metadata/reference_metadata.tsv

    # Put sendsketch's output into a single folder for organization
    mkdir ${prefix}_inputs/sendsketch
    cp -r ${sendsketch} ${prefix}_inputs/sendsketch/

    # Put all of the statistics for NCBI references considered into a single folder
    if [ ! -z "${ncbi_ref_meta}" ]; then
        mkdir -p ${prefix}_inputs/reference_data/considered
        cp -r ${ncbi_ref_meta} ${prefix}_inputs/reference_data/considered/
    fi

    # Put the metadata for references selected for each sample into a single folder
    if [ ! -z "${selected_refs}" ]; then
        mkdir -p ${prefix}_inputs/reference_data/downloaded
        cp -r ${selected_refs} ${prefix}_inputs/reference_data/downloaded/
    fi

    # Add estimated ANI matrix from sourmash
    mkdir -p ${prefix}_inputs/sketch_comparisons/ani_matricies
    cp ${ani_matrix} ${prefix}_inputs/sketch_comparisons/ani_matricies/

    # Add metadata for references assined for variant calling
    if [ ! -z "${mapping_ref}" ]; then
        mkdir -p ${prefix}_inputs/reference_data/selected
        cp ${mapping_ref} ${prefix}_inputs/reference_data/selected/${prefix}_mapping_references.tsv
    fi

    # Add SNP alignment from variant calling
    if [ ! -z "${snp_aligns}" ]; then
        mkdir -p ${prefix}_inputs/variants
        cp -r ${snp_aligns} ${prefix}_inputs/variants
    fi

    # Add SNP phylogenies from variant calling
    if [ ! -z "${snp_phylos}" ]; then
        mkdir -p ${prefix}_inputs/trees/snp
        cp -r ${snp_phylos} ${prefix}_inputs/trees/snp/
    fi

    # Add seleted references for the core gene phylogeny
    if [ ! -z "${core_phylo_refs}" ]; then
        mkdir -p ${prefix}_inputs/reference_data/selected/
        cp -r ${core_phylo_refs} ${prefix}_inputs/reference_data/selected/${prefix}_core_references.tsv
    fi

    # Add POCP estimate from the core genome phylogeny
    if [ ! -z "${pocp}" ]; then
        mkdir -p ${prefix}_inputs/pocp
        cp -r ${pocp} ${prefix}_inputs/pocp/${prefix}_pocp.tsv
    fi

    # Add core genome phylogenies
    if [ ! -z "${core_phylos}" ]; then
        mkdir -p ${prefix}_inputs/trees/core
        cp -r ${core_phylos} ${prefix}_inputs/trees/core/
    fi

    # Add seleted references for the busco phylogeny
    if [ ! -z "${busco_refs}" ]; then
        mkdir -p ${prefix}_inputs/reference_data/selected/
        cp -r ${busco_refs} ${prefix}_inputs/reference_data/selected/${prefix}_busco_references.tsv
    fi

    # Add busco phylogeny
    if [ ! -z "${busco_phylos}" ]; then
        mkdir -p ${prefix}_inputs/trees/busco
        cp -r ${busco_phylos} ${prefix}_inputs/trees/busco/
    fi

    # Put multiqc's output into a single folder for organization
    cp -r ${multiqc} ${prefix}_inputs/multiqc

    # Add pipeline status messages
    mkdir -p ${prefix}_inputs/pipeline_info
    if [ -z "${messages}" ]; then
        printf '"report_id"\t"sample_id"\t"reference_id"\t"workflow"\t"level"\t"message"\\n' > ${prefix}_inputs/pipeline_info/messages.tsv
    else
        cp ${messages} ${prefix}_inputs/pipeline_info/messages.tsv
    fi

    # Add versions of software used
    mkdir -p ${prefix}_inputs/pipeline_info
    cp ${versions} ${prefix}_inputs/pipeline_info/version_info.yml

    # Record run metadata
    mkdir -p ${prefix}_inputs/pipeline_info
    cat <<-END_INFO > ${prefix}_inputs/pipeline_info/pathogensurveillance_run_info.yml
    group_id: ${group_meta.id}
    container_engine: ${workflow.containerEngine}
    profile: ${workflow.profile}
    nextflow_version: ${nextflow.version}
    pipeine_verison: ${workflow.manifest.version}
    END_INFO
    """
}
