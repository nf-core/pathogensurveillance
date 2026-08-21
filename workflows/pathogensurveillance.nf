/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { PREPARE_INPUT               } from '../subworkflows/local/prepare_input'
include { CORE_GENOME_PHYLOGENY       } from '../subworkflows/local/core_genome_phylogeny'
include { VARIANT_ANALYSIS            } from '../subworkflows/local/variant_analysis'
include { SKETCH_COMPARISON           } from '../subworkflows/local/sketch_comparison'
include { GENOME_ASSEMBLY             } from '../subworkflows/local/genome_assembly'
include { BUSCO_PHYLOGENY             } from '../subworkflows/local/busco_phylogeny'
include { INITIAL_QC_CHECKS           } from '../subworkflows/local/initial_qc_checks'
include { MAIN_REPORT                 } from '../modules/local/main_report'
include { DOWNLOAD_ASSEMBLIES         } from '../modules/local/download_assemblies'
include { PREPARE_REPORT_INPUT        } from '../modules/local/prepare_report_input'
include { softwareVersionsToYAML      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { paramsSummaryMap            } from 'plugin/nf-schema'
include { paramsSummaryMultiqc        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText      } from '../subworkflows/local/utils_nfcore_pathogensurveillance_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PATHOGENSURVEILLANCE {

    take:
    sample_data_tsv
    reference_data_tsv

    main:

    // Write output format file for PathoSurveilR parsing
    file("$projectDir/assets/.pathogensurveillance_output.json").copyTo("${params.outdir}/.pathogensurveillance_output.json")

    // Initalize channel to accumulate information about software versions used
    versions = channel.empty()
    messages = channel.empty()

    // Read in samplesheet, validate and stage input files
    PREPARE_INPUT ( sample_data_tsv, reference_data_tsv )
    versions = versions.mix(PREPARE_INPUT.out.versions)
    messages = messages.mix(PREPARE_INPUT.out.messages)

    // Assemble and annotate genomes
    GENOME_ASSEMBLY (
        PREPARE_INPUT.out.sample_data
    )
    versions = versions.mix(GENOME_ASSEMBLY.out.versions)
    messages = messages.mix(GENOME_ASSEMBLY.out.messages)

    // Initial quick analysis of sequences and references based on sketchs
    SKETCH_COMPARISON (
        PREPARE_INPUT.out.sample_data,
        GENOME_ASSEMBLY.out.scaffolds
    )
    versions = versions.mix(SKETCH_COMPARISON.out.versions)
    messages = messages.mix(SKETCH_COMPARISON.out.messages)

    // Initial quality control of reads
    INITIAL_QC_CHECKS ( PREPARE_INPUT.out.sample_data )
    versions = versions.mix(INITIAL_QC_CHECKS.out.versions)
    messages = messages.mix(INITIAL_QC_CHECKS.out.messages)

    // Call variants and create SNP-tree and minimum spanning nextwork
    VARIANT_ANALYSIS (
        PREPARE_INPUT.out.sample_data,
        SKETCH_COMPARISON.out.ani_matrix
    )
    versions = versions.mix(VARIANT_ANALYSIS.out.versions)
    messages = messages.mix(VARIANT_ANALYSIS.out.messages)

    // Create core gene phylogeny for bacterial samples
    if (!params.skip_core_phylogeny) {
        CORE_GENOME_PHYLOGENY (
            PREPARE_INPUT.out.sample_data,
            SKETCH_COMPARISON.out.ani_matrix,
            GENOME_ASSEMBLY.out.scaffolds
        )
        versions = versions.mix(CORE_GENOME_PHYLOGENY.out.versions)
        messages  = messages.mix(CORE_GENOME_PHYLOGENY.out.messages)
        core_selected_refs = CORE_GENOME_PHYLOGENY.out.selected_refs
        core_pocp = CORE_GENOME_PHYLOGENY.out.pocp
        core_phylogeny = CORE_GENOME_PHYLOGENY.out.phylogeny
    } else {
        core_selected_refs = channel.empty()
        core_pocp = channel.empty()
        core_phylogeny = channel.empty()
    }

    // Read2tree BUSCO phylogeny for eukaryotes
    BUSCO_PHYLOGENY (
        PREPARE_INPUT.out.sample_data,
        SKETCH_COMPARISON.out.ani_matrix,
        GENOME_ASSEMBLY.out.scaffolds
    )
    versions = versions.mix(BUSCO_PHYLOGENY.out.versions)
    messages = messages.mix(BUSCO_PHYLOGENY.out.messages)

    // Collate and save software versions
    def topic_versions_all = channel.topic("versions")
        .distinct()
    def topic_versions_file = topic_versions_all.filter { entry -> entry instanceof Path }
    def topic_versions_tuple = topic_versions_all.filter { entry -> !(entry instanceof Path) }

    def topic_versions_string = topic_versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(versions.mix(topic_versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'version_info.yml',
            sort: true,
            newLine: true
         ).set { collated_versions }

    // MultiQC
    multiqc_config          = channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    multiqc_custom_config   = params.multiqc_config ? channel.fromPath( params.multiqc_config, checkIfExists: true ) : channel.empty()
    multiqc_logo            = params.multiqc_logo   ? channel.fromPath( params.multiqc_logo, checkIfExists: true ) : channel.empty()
    multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    // Note: the belowsection was from a template update that has not been merged into this logic yet
    methods_description     = channel.value(methodsDescriptionText(multiqc_custom_methods_description))
    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    // End note section -------------------

    fastqc_results = PREPARE_INPUT.out.sample_data
        .map{ sample_meta -> [[id: sample_meta.sample_id], [id: sample_meta.report_group_ids]] }
        .combine(INITIAL_QC_CHECKS.out.fastqc_zip, by: 0)
        .map{ sample_meta, report_meta, fastqc -> [report_meta, fastqc] }
        .unique()
        .groupTuple(sort: 'hash')
    fastp_results = PREPARE_INPUT.out.sample_data
        .map{ sample_meta -> [[id: sample_meta.sample_id], [id: sample_meta.report_group_ids]] }
        .combine(GENOME_ASSEMBLY.out.fastp_json, by: 0)
        .map{ sample_meta, report_meta, fastp_json -> [report_meta, fastp_json] }
        .unique()
        .groupTuple(sort: 'hash')
    nanoplot_results = PREPARE_INPUT.out.sample_data
        .map{ sample_meta -> [[id: sample_meta.sample_id], [id: sample_meta.report_group_ids]] }
        .combine(INITIAL_QC_CHECKS.out.nanoplot_txt, by: 0)
        .map{ sample_meta, report_meta, nanoplot_txt -> [report_meta, nanoplot_txt] }
        .unique()
        .groupTuple(sort: 'hash')
    quast_results = PREPARE_INPUT.out.sample_data
        .map{ sample_meta -> [[id: sample_meta.sample_id], [id: sample_meta.report_group_ids]] }
        .combine(GENOME_ASSEMBLY.out.quast, by: 0)
        .map{ sample_meta, report_meta, quast -> [report_meta, quast] }
        .unique()
        .groupTuple(sort: 'hash')
    multiqc_files = PREPARE_INPUT.out.sample_data
        .map{ sample_meta -> [[id: sample_meta.report_group_ids]] }
        .unique()
        .combine(collated_versions)
        .join(fastqc_results, remainder: true)
        .join(fastp_results, remainder: true)
        .join(nanoplot_results, remainder: true)
        .join(quast_results, remainder: true)
        .map {report_meta, my_versions, fastqc, fastp, nanoplot, quast ->
            def files = (fastqc ?: []) + (fastp ?: []) + (nanoplot ?: []) + (quast ?: []) + ([my_versions])
            [report_meta, files.flatten()]
        }
    def multiqc_all = multiqc_files
        .combine(multiqc_config.collect(sort: true).ifEmpty([]))
        .combine(multiqc_custom_config.collect(sort: true).ifEmpty([]))
        .combine(multiqc_logo.collect(sort: true).ifEmpty([]))
        .combine(channel.value([]))
        .combine(channel.value([]))
        .map { tuple ->
            def report_meta = tuple[0]
            def files = tuple[1]
            def config = tuple.size() > 2 ? tuple[2] : []
            def custom_config = tuple.size() > 3 ? tuple[3] : []
            def logo = tuple.size() > 4 ? tuple[4] : []
            def replace = tuple.size() > 5 ? tuple[5] : []
            def samples = tuple.size() > 6 ? tuple[6] : []
            def all_configs = (config ? [config] : []) + (custom_config ? [custom_config] : [])
            [report_meta, files, all_configs.flatten(), logo, replace, samples]
        }

    MULTIQC ( multiqc_all )
    versions = versions.mix(MULTIQC.out.versions)

    // Gather sample data for each report
    sample_data_tsvs = PREPARE_INPUT.out.sample_data
        .map{ sample_meta ->
            [[id: sample_meta.report_group_ids], sample_meta.findAll { entry -> entry.key != 'paths' && entry.key != 'ref_metas' && entry.key != 'ref_ids' }]
        }
        .unique()
        .collectFile(keepHeader: true, skip: 1) { report_meta, sample_meta ->
            [ "${report_meta.id}_sample_data.tsv", sample_meta.keySet().collect{ key -> '"' + key + '"'}.join('\t') + "\n" + sample_meta.values().collect{ value -> '"' + (value ?: '') + '"'}.join('\t') + "\n" ]
        }
        .map { file ->[[id: file.getSimpleName().replace('_sample_data', '')], file]}

    // Gather reference data for each report
    reference_data_tsvs = PREPARE_INPUT.out.sample_data
        .map { sample_meta ->
            [[id: sample_meta.report_group_ids], sample_meta.ref_metas]
        }
        .transpose(by: 1)
        .map { report_meta, ref_meta ->
            [report_meta, ref_meta.findAll { entry -> entry.key != 'ref_path' && entry.key != 'gff' }]
        }
        .unique()
        .collectFile(keepHeader: true, skip: 1) { report_meta, ref_meta ->
            [ "${report_meta.id}_reference_data.tsv", ref_meta.keySet().collect{ key -> '"' + key + '"'}.join('\t') + "\n" + ref_meta.values().collect{ value -> '"' + (value ?: '') + '"'}.join('\t') + "\n" ]
        }
        .map { file ->[[id: file.getSimpleName().replace('_reference_data', '')], file]}

    // Gather sendsketch signatures and taxa found
    sendsketch_files = PREPARE_INPUT.out.sample_data
        .map{ sample_meta -> [[id: sample_meta.sample_id], [id: sample_meta.report_group_ids]] }
        .combine(PREPARE_INPUT.out.sendsketch, by: 0)
        .map{ sample_meta, report_meta, sendsketch -> [report_meta, sendsketch] }
        .unique()
    sendsketch_taxa = PREPARE_INPUT.out.sample_data
        .map{ sample_meta -> [[id: sample_meta.sample_id], [id: sample_meta.report_group_ids]] }
        .combine(PREPARE_INPUT.out.taxa_found, by: 0)
        .map{ sample_meta, report_meta, taxa_found -> [report_meta, taxa_found] }
        .unique()
    sendsketch_hits = sendsketch_files
        .mix(sendsketch_taxa)
        .groupTuple(by: 0, sort: 'hash')

    // Gather NCBI reference metadata for all references considered
    ncbi_ref_meta = PREPARE_INPUT.out.sample_data
        .map{ sample_meta -> [[id: sample_meta.sample_id], [id: sample_meta.report_group_ids]] }
        .combine(PREPARE_INPUT.out.family_stats_per_sample, by: 0)
        .groupTuple(by: 1, sort: 'hash')
        .map { sample_meta, report_meta, family_stats ->
            [report_meta, family_stats.flatten().unique()]
        }

    // Gather selected reference metadata
    selected_ref_meta = PREPARE_INPUT.out.sample_data
        .map{ sample_meta -> [[id: sample_meta.sample_id], [id: sample_meta.report_group_ids]] }
        .combine(PREPARE_INPUT.out.selected_ref_meta, by:0)
        .map{ sample_meta, report_meta, ref_meta_file ->
            [report_meta, ref_meta_file] }
        .unique()
        .groupTuple(sort: 'hash')
        .map { report_meta, ref_meta_files ->
            [report_meta, ref_meta_files.findAll{ file -> file != null }]
        }

    // Gather SNP alignments from the variant analysis
    snp_align = VARIANT_ANALYSIS.out.snp_align
        .map { report_meta, ref_meta, fasta -> [report_meta, fasta] }
        .groupTuple(sort: 'hash')

    // Gather phylogenies from the variant analysis
    snp_phylogeny = VARIANT_ANALYSIS.out.phylogeny
        .map { report_meta, ref_meta, tree -> [report_meta, tree] }
        .groupTuple(sort: 'hash')

    // Gather status messages for each group
    group_messages = messages
        .unique()
        .collectFile(keepHeader: true, skip: 1) { sample_meta, report_meta, ref_meta, workflow, level, message ->
            [ "${report_meta.id}.tsv", "\"report_id\"\t\"sample_id\"\t\"reference_id\"\t\"workflow\"\t\"level\"\t\"message\"\n\"${report_meta.id}\"\t\"${sample_meta ? sample_meta.id : ''}\"\t\"${ref_meta ? ref_meta.id : ''}\"\t\"${workflow}\"\t\"${level}\"\t\"${message}\"\n" ]
        }
        .map { file ->[[id: file.getSimpleName()], file]}
        .ifEmpty([])

    // Combine components into a single channel for the main report_meta
    report_inputs = sample_data_tsvs
        .join(reference_data_tsvs, remainder: true)
        .join(sendsketch_hits, remainder: true)
        .join(ncbi_ref_meta, remainder: true)
        .join(selected_ref_meta, remainder: true)
        .join(VARIANT_ANALYSIS.out.ani_matrix, remainder: true)
        .join(VARIANT_ANALYSIS.out.mapping_ref, remainder: true)
        .join(snp_align, remainder: true)
        .join(snp_phylogeny, remainder: true)
        .join(core_selected_refs, remainder: true)
        .join(core_pocp, remainder: true)
        .join(core_phylogeny, remainder: true)
        .join(BUSCO_PHYLOGENY.out.selected_refs, remainder: true)
        .join(BUSCO_PHYLOGENY.out.tree, remainder: true)
        .join(MULTIQC.out.report, remainder: true)
        .join(group_messages, remainder: true)
        .filter{ item -> item[0] != null }
        .map{ item -> item.size() == 16 ? item + [null] : item }
        .filter{ item -> item.size() == 17 }
        .map{ item -> item.collect{ element -> element ?: [] } }
        .combine(collated_versions)

    PREPARE_REPORT_INPUT (
        report_inputs,
        channel.fromPath("${projectDir}/assets/.pathogensurveillance_output.json", checkIfExists: true).first()
    )

    MAIN_REPORT (
        PREPARE_REPORT_INPUT.out.report_input,
        channel.fromPath("${projectDir}/assets/main_report", checkIfExists: true).first()
    )

    // Collate and save messages
    messages
        .unique()
        .map  { sample_meta, report_meta, ref_meta, workflow, level, message ->
            "\"report_id\"\t\"sample_id\"\t\"reference_id\"\t\"workflow\"\t\"level\"\t\"message\"\n\"${report_meta.id}\"\t\"${sample_meta ? sample_meta.id : ''}\"\t\"${ref_meta ? ref_meta.id : ''}\"\t\"${workflow}\"\t\"${level}\"\t\"${message}\"\n"
        }
        .ifEmpty("\"report_id\"\t\"sample_id\"\t\"reference_id\"\t\"workflow\"\t\"level\"\t\"message\"\n")
        .collectFile(
            keepHeader: true,
            skip: 1,
            storeDir: "${params.outdir}/pipeline_info",
            name: "messages.tsv",
            sort: true
        )

    // Save pipeline execution paramters
    channel.value(
        """
        command_line: ${workflow.commandLine}
        commit_id: ${workflow.commitId}
        container_engine: ${workflow.containerEngine}
        profile: ${workflow.profile}
        revision: ${workflow.revision}
        run_name: ${workflow.runName}
        session_id: ${workflow.sessionId}
        start_time: ${workflow.start}
        nextflow_version: ${nextflow.version}
        pipeline_version: ${workflow.manifest.version}
        """.stripIndent().trim()
    )
    .collectFile(storeDir: "${params.outdir}/pipeline_info", name: "pathogensurveillance_run_info.yml")

    // Gather sample data for each report
    PREPARE_INPUT.out.sample_data
        .map{ sample_meta ->
            sample_meta.findAll { entry -> entry.key != 'paths' && entry.key != 'ref_metas' && entry.key != 'ref_ids' }
        }
        .unique()
        .collectFile(
            keepHeader: true,
            skip: 1,
            storeDir: "${params.outdir}/metadata",
            name: "sample_metadata.tsv"
        ) { sample_meta ->
            sample_meta.keySet().collect{ key -> '"' + key + '"'}.join('\t') + "\n" + sample_meta.values().collect{ value -> '"' + (value ?: '') + '"'}.join('\t') + "\n"
        }

    // Gather reference data for each report
    reference_data_tsvs = PREPARE_INPUT.out.sample_data
        .map { sample_meta ->
            [sample_meta.ref_metas]
        }
        .transpose(by: 0)
        .map { ref_meta ->
            ref_meta[0].findAll { entry -> entry.key != 'ref_path' && entry.key != 'gff' }
        }
        .unique()
        .collectFile(
            keepHeader: true,
            skip: 1,
            storeDir: "${params.outdir}/metadata",
            name: "reference_metadata.tsv"
        ) { ref_meta ->
            ref_meta.keySet().collect{ key -> '"' + key + '"'}.join('\t') + "\n" + ref_meta.values().collect{ value -> '"' + (value ?: '') + '"'}.join('\t') + "\n"
        }

    emit:
    multiqc_report = MULTIQC.out.report
    versions       = versions                 // channel: [ path(versions.yml) ]
}
