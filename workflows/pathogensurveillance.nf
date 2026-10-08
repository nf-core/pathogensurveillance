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

// Built-in template directory aliases. The default directory was renamed to pathsurveil_report,
// but "report" is what users have always written in report_data, so it keeps working. The
// canonical directory name is also what names the published report, so the alias has to be
// applied before the value reaches MAIN_REPORT, not merely when resolving the path. Kept in
// step with the alias map in bin/check_samplesheet.R.
def templateAliases() {
    return ['report': 'pathsurveil_report']
}

// Canonical form of a report_data template spec: aliases are replaced, absolute paths are left
// alone (there is nothing to alias, and their basename becomes the name of the report).
def canonicalTemplateSpec(String spec) {
    def value = spec?.toString()?.trim() ?: ''
    if (value.startsWith('/')) {
        return value
    }
    def aliases = templateAliases()
    return aliases.containsKey(value) ? aliases[value] : value
}

def resolveTemplateDir(String spec) {
    // Kept in step with template_status() in bin/check_samplesheet.R: a bare name resolves under
    // assets/report_templates/, an absolute path is used as-is, and anything else is rejected.
    // Relative paths must never be resolved here. file() would anchor them to the launch
    // directory rather than the project root, while the R validator resolves against the task
    // work dir, so the two would disagree on the same input.
    def value = canonicalTemplateSpec(spec)
    if (value.startsWith('~')) {
        error("report_data template '${value}' starts with '~', which is not expanded. Use an absolute path (e.g. /data/templates/tpl) or a name resolved under assets/report_templates/.")
    }
    if (value.startsWith('/')) {
        return file(value)
    }
    if (value.contains('/') || value.startsWith('.')) {
        error("report_data template '${value}' is a relative path. Use an absolute path (e.g. /data/templates/tpl) or a name resolved under assets/report_templates/.")
    }
    return file("${projectDir}/assets/report_templates/${value}")
}

// Normalise a field read from a metadata table: trim and drop surrounding quotes so that
// report_group_ids/template line up with the unquoted values coming from the samplesheet.
def cleanField(value) {
    def text = value?.toString() ?: ''
    text = text.replaceAll('^\s*["\']', '').replaceAll('["\']\s*$', '')
    text.trim()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PATHOGENSURVEILLANCE {

    take:
    sample_data_tsv
    reference_data_tsv
    report_data_tsv

    main:

    // Write output format file for PathoSurveilR parsing
    file("$projectDir/assets/.pathogensurveillance_output.json").copyTo("${params.outdir}/.pathogensurveillance_output.json")

    // Initialize channel to accumulate warning messages
    messages = channel.empty()

    // Read in samplesheet, validate and stage input files
    PREPARE_INPUT ( sample_data_tsv, reference_data_tsv, report_data_tsv )
    messages = messages.mix(PREPARE_INPUT.out.messages)

    // Assemble and annotate genomes
    GENOME_ASSEMBLY (
        PREPARE_INPUT.out.sample_data
    )
    messages = messages.mix(GENOME_ASSEMBLY.out.messages)

    // Initial quick analysis of sequences and references based on sketchs
    SKETCH_COMPARISON (
        PREPARE_INPUT.out.sample_data,
        GENOME_ASSEMBLY.out.scaffolds
    )
    messages = messages.mix(SKETCH_COMPARISON.out.messages)

    // Initial quality control of reads
    INITIAL_QC_CHECKS ( PREPARE_INPUT.out.sample_data )
    messages = messages.mix(INITIAL_QC_CHECKS.out.messages)

    // Call variants and create SNP-tree and minimum spanning nextwork
    VARIANT_ANALYSIS (
        PREPARE_INPUT.out.sample_data,
        SKETCH_COMPARISON.out.pairwise_csv
    )
    messages = messages.mix(VARIANT_ANALYSIS.out.messages)

    // Create core gene phylogeny for bacterial samples
    if (!params.skip_core_phylogeny) {
        CORE_GENOME_PHYLOGENY (
            PREPARE_INPUT.out.sample_data,
            SKETCH_COMPARISON.out.pairwise_csv,
            GENOME_ASSEMBLY.out.scaffolds
        )
        messages  = messages.mix(CORE_GENOME_PHYLOGENY.out.messages)
        core_selected_refs = CORE_GENOME_PHYLOGENY.out.selected_refs
        core_pocp = CORE_GENOME_PHYLOGENY.out.pocp
        core_phylogeny = CORE_GENOME_PHYLOGENY.out.phylogeny
        core_gene_counts = CORE_GENOME_PHYLOGENY.out.gene_count
    } else {
        core_selected_refs = channel.empty()
        core_pocp = channel.empty()
        core_phylogeny = channel.empty()
        core_gene_counts = channel.empty()
    }

    // Read2tree BUSCO phylogeny for eukaryotes
    BUSCO_PHYLOGENY (
        PREPARE_INPUT.out.sample_data,
        SKETCH_COMPARISON.out.pairwise_csv,
        GENOME_ASSEMBLY.out.scaffolds
    )
    messages = messages.mix(BUSCO_PHYLOGENY.out.messages)
    busco_gene_counts = BUSCO_PHYLOGENY.out.gene_count

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

    softwareVersionsToYAML(topic_versions_file)
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

    // Gather gene counts for report
    core_gene_counts_grouped = core_gene_counts
        .map { report_meta, tsv -> [report_meta, tsv] }
        .groupTuple(sort: 'hash')
    busco_gene_counts_grouped = busco_gene_counts
        .map { report_meta, tsv -> [report_meta, tsv] }
        .groupTuple(sort: 'hash')
    gene_counts = core_gene_counts_grouped
        .mix(busco_gene_counts_grouped)
        .groupTuple(sort: 'hash')
        .map { report_meta, tsvs -> [report_meta, tsvs.flatten().unique()] }

    // Resolve report templates per group (report_data optional, default "report")
    // PREPARE_INPUT.out.report_data is a TSV file (header report_group_ids,template)
    // Items are [report_meta, [template, is_explicit]]; is_explicit marks a template that came
    // from report_data rather than the built-in default. The pair is nested in a list because
    // groupTuple flattens the non-key elements of each item.
    report_template_specs = PREPARE_INPUT.out.report_data
        .splitCsv(header: true, sep: '\t', quote: '"')
        .map { row -> [cleanField(row.report_group_ids), cleanField(row.template)] }
        .flatMap { entry ->
            def report_ids = entry[0]
            def tmpl = entry[1]
            def groups = report_ids.toString().split(';').collect{ it.trim() }.findAll{ it }
            def templates = tmpl.toString().split(';').collect{ it.trim() }.findAll{ it }
            if (groups.isEmpty() || templates.isEmpty()) return []
            groups.collectMany{ g -> templates.collect{ t -> [[id: g], [t, true]] } }
        }

    // Every group gets the default "report" template. Defaults and report_data entries are
    // mixed into ONE channel and grouped by key, so the decision is made per group without
    // join()/combine() -- both of which silently drop every group when report_data is absent.
    default_per_group = PREPARE_INPUT.out.sample_data
        .map { sample_meta -> [[id: sample_meta.report_group_ids]] }
        .unique()
        .map { report_meta -> [report_meta, ['report', false]] }

    template_per_group = report_template_specs
        .mix(default_per_group)
        .groupTuple(by: 0, sort: 'hash')
        .flatMap { entry ->
            def report_meta = entry[0]
            def specs = entry[1]
            // An explicit report_data template shadows the default for that group, so a group
            // listing only a non-report template renders only that template, while a group
            // listing "report" plus others renders all of them (deduplicated).
            def explicit = specs.findAll { spec -> spec[1] }.collect { spec -> spec[0] }.unique()
            def templates = explicit ?: ['report']
            templates.collect { tmpl -> [report_meta, tmpl] }
        }
    template_dirs = template_per_group
        .map { pair ->
            // Canonicalise before this becomes group_meta.template: MAIN_REPORT derives the
            // published filename from it, so an unaliased "report" would name the file
            // all_report.html instead of all_pathsurveil_report.html.
            def canon = canonicalTemplateSpec(pair[1])
            [pair[0], canon, resolveTemplateDir(canon)]
        }

    // Combine components into a single channel for the main report_meta
    report_inputs = sample_data_tsvs
        .join(reference_data_tsvs, remainder: true)
        .join(sendsketch_hits, remainder: true)
        .join(ncbi_ref_meta, remainder: true)
        .join(selected_ref_meta, remainder: true)
        .join(SKETCH_COMPARISON.out.pairwise_csv, remainder: true)
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
        .join(gene_counts, remainder: true)
        .filter{ item -> item[0] != null }
        .map{ item -> item.size() == 17 ? item + [null] : item }
        .filter{ item -> item.size() == 18 }
        .map{ item -> item.collect{ element -> element ?: [] } }
        .combine(collated_versions)

    PREPARE_REPORT_INPUT (
        report_inputs,
        channel.fromPath("${projectDir}/assets/.pathogensurveillance_output.json", checkIfExists: true).first()
    )

    // Combine per-report inputs with per-template dirs for rendering
    // Each report group may have multiple templates (e.g. report;dashboard) -> one render per (group,template)
    // combine(by: 0) returns [key, report_input, template, template_dir] -- the key element of each
    // source becomes element 0 and the remaining elements of each source follow in source order.
    def main_report_inputs = PREPARE_REPORT_INPUT.out.report_input
        .combine(template_dirs, by: 0)
        .map { entry -> [[id: entry[0].id, template: entry[2]], entry[1], entry[3]] }

    MAIN_REPORT(
        main_report_inputs
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
}
