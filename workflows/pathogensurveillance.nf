/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC                     } from '../modules/local/multiqc/main'
include { PREPARE_INPUT               } from '../subworkflows/local/prepare_input'
include { CORE_GENOME_PHYLOGENY       } from '../subworkflows/local/core_genome_phylogeny'
include { VARIANT_ANALYSIS            } from '../subworkflows/local/variant_analysis'
include { SKETCH_COMPARISON           } from '../subworkflows/local/sketch_comparison'
include { GENOME_ASSEMBLY             } from '../subworkflows/local/genome_assembly'
include { BUSCO_PHYLOGENY             } from '../subworkflows/local/busco_phylogeny'
include { INITIAL_QC_CHECKS           } from '../subworkflows/local/initial_qc_checks'
include { MAIN_REPORT                 } from '../modules/local/main_report'
include { RECORD_MESSAGES             } from '../modules/local/record_messages'
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


    // Initalize channel to accumulate information about software versions used
    versions = Channel.empty()
    messages = Channel.empty()

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
    CORE_GENOME_PHYLOGENY (
        PREPARE_INPUT.out.sample_data,
        SKETCH_COMPARISON.out.ani_matrix,
        GENOME_ASSEMBLY.out.scaffolds
    )
    versions = versions.mix(CORE_GENOME_PHYLOGENY.out.versions)
    messages  = messages.mix(CORE_GENOME_PHYLOGENY.out.messages)

    // Read2tree BUSCO phylogeny for eukaryotes
    BUSCO_PHYLOGENY (
        PREPARE_INPUT.out.sample_data,
        SKETCH_COMPARISON.out.ani_matrix,
        GENOME_ASSEMBLY.out.scaffolds
    )
    versions = versions.mix(BUSCO_PHYLOGENY.out.versions)
    messages = messages.mix(BUSCO_PHYLOGENY.out.messages)

    // Collate and save software versions
    softwareVersionsToYAML(versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { collated_versions }

    // MultiQC
    multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
    multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
    multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    fastqc_results = PREPARE_INPUT.out.sample_data
        .map{ [[id: it.sample_id], [id: it.report_group_ids]] }
        .combine(INITIAL_QC_CHECKS.out.fastqc_zip, by: 0)
        .map{ sample_meta, report_meta, fastqc -> [report_meta, fastqc] }
        .unique()
        .groupTuple(sort: 'hash')
    nanoplot_results = PREPARE_INPUT.out.sample_data
        .map{ [[id: it.sample_id], [id: it.report_group_ids]] }
        .combine(INITIAL_QC_CHECKS.out.nanoplot_txt, by: 0)
        .map{ sample_meta, report_meta, nanoplot_txt -> [report_meta, nanoplot_txt] }
        .unique()
        .groupTuple(sort: 'hash')
    quast_results = PREPARE_INPUT.out.sample_data
        .map{ [[id: it.sample_id], [id: it.report_group_ids]] }
        .combine(GENOME_ASSEMBLY.out.quast, by: 0)
        .map{ sample_meta, report_meta, quast -> [report_meta, quast] }
        .unique()
        .groupTuple(sort: 'hash')
    multiqc_files = fastqc_results
        .join(nanoplot_results, remainder: true)
        .join(quast_results, remainder: true)
        .combine(collated_versions)
        .map { report_meta, fastqc, nanoplot, quast, versions ->
            files = fastqc ?: [] + nanoplot ?: [] + quast ?: [] + [versions]
            [report_meta, files.flatten()]
        }
    MULTIQC (
        multiqc_files,
        multiqc_config.collect(sort: true).ifEmpty([]),
        multiqc_custom_config.collect(sort: true).ifEmpty([]),
        multiqc_logo.collect(sort: true).ifEmpty([]),
        [],
        []
    )
    versions = versions.mix(MULTIQC.out.versions)

    // Gather sample data for each report
    sample_data_tsvs = PREPARE_INPUT.out.sample_data
        .map{ sample_meta ->
            [[id: sample_meta.report_group_ids], sample_meta.findAll {it.key != 'paths' && it.key != 'ref_metas' && it.key != 'ref_ids'}]
        }
        .unique()
        .collectFile(keepHeader: true, skip: 1) { report_meta, sample_meta ->
            [ "${report_meta.id}_sample_data.tsv", sample_meta.keySet().collect{'"' + it + '"'}.join('\t') + "\n" + sample_meta.values().collect{'"' + it + '"'}.join('\t') + "\n" ]
        }
        .map {[[id: it.getSimpleName().replace('_sample_data', '')], it]}

    // Gather reference data for each report
    reference_data_tsvs = PREPARE_INPUT.out.sample_data
        .map { sample_meta ->
            [[id: sample_meta.report_group_ids], sample_meta.ref_metas]
        }
        .transpose(by: 1)
        .map { report_meta, ref_meta ->
            [report_meta, ref_meta.findAll {it.key != 'ref_path' && it.key != 'gff'}]
        }
        .unique()
        .collectFile(keepHeader: true, skip: 1) { report_meta, ref_meta ->
            [ "${report_meta.id}_reference_data.tsv", ref_meta.keySet().collect{'"' + it + '"'}.join('\t') + "\n" + ref_meta.values().collect{'"' + it + '"'}.join('\t') + "\n" ]
        }
        .map {[[id: it.getSimpleName().replace('_reference_data', '')], it]}

    // Gather sendsketch signatures
    sendsketch_hits = PREPARE_INPUT.out.sample_data
        .map{ [[id: it.sample_id], [id: it.report_group_ids]] }
        .combine(PREPARE_INPUT.out.sendsketch, by: 0)
        .map{ sample_meta, report_meta, sendsketch -> [report_meta, sendsketch] }
        .unique()
        .groupTuple(sort: 'hash')

    // Gather NCBI reference metadata for all references considered
    ncbi_ref_meta = PREPARE_INPUT.out.sample_data
        .map{ [[id: it.sample_id], [id: it.report_group_ids]] }
        .combine(PREPARE_INPUT.out.family_stats_per_sample, by: 0)
        .groupTuple(by: 1, sort: 'hash')
        .map { sample_meta, report_meta, family_stats ->
            [report_meta, family_stats.flatten().unique()]
        }

    // Gather selected reference metadata
    selected_ref_meta = PREPARE_INPUT.out.sample_data
        .map{ [[id: it.sample_id], [id: it.report_group_ids]] }
        .combine(PREPARE_INPUT.out.selected_ref_meta, by:0)
        .map{ sample_meta, report_meta, ref_meta_file ->
            [report_meta, ref_meta_file] }
        .unique()
        .groupTuple(sort: 'hash')
        .map { report_meta, ref_meta_files ->
            [report_meta, ref_meta_files.findAll{it != null}]
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
            [ "${report_meta.id}.tsv", "\"sample_id\"\t\"reference_id\"\t\"workflow\"\t\"level\"\t\"message\"\n\"${sample_meta ? sample_meta.id : 'NA'}\"\t\"${ref_meta ? ref_meta.id : 'NA'}\"\t\"${workflow}\"\t\"${level}\"\t\"${message}\"\n" ]
        }
        .map {[[id: it.getSimpleName()], it]}
        .ifEmpty([])

    // Combine components into a single channel for the main report_meta
    report_inputs = sample_data_tsvs
        .join(reference_data_tsvs, remainder: true)
        .join(sendsketch_hits, remainder: true)
        .join(ncbi_ref_meta, remainder: true)
        .join(selected_ref_meta, remainder: true)
        .join(SKETCH_COMPARISON.out.ani_matrix, remainder: true)
        .join(VARIANT_ANALYSIS.out.mapping_ref, remainder: true)
        .join(snp_align, remainder: true)
        .join(snp_phylogeny, remainder: true)
        .join(CORE_GENOME_PHYLOGENY.out.selected_refs, remainder: true)
        .join(CORE_GENOME_PHYLOGENY.out.pocp, remainder: true)
        .join(CORE_GENOME_PHYLOGENY.out.phylogeny, remainder: true)
        .join(BUSCO_PHYLOGENY.out.selected_refs, remainder: true)
        .join(BUSCO_PHYLOGENY.out.tree, remainder: true)
        .join(MULTIQC.out.outdir, remainder: true)
        .join(group_messages, remainder: true)
        .filter{it[0] != null} // remove extra item if messages is empty
        .map{ it.size() == 16 ? it + [null] : it } // adds placeholder if messages is empty
        .map{ it.collect{ it ?: [] } } //replace nulls with empty lists

    PREPARE_REPORT_INPUT (
        report_inputs,
        collated_versions.first() // .first converts it to a value channel so it can be reused for multiple reports.
    )

    MAIN_REPORT (
        PREPARE_REPORT_INPUT.out.report_input,
        Channel.fromPath("${projectDir}/assets/main_report", checkIfExists: true).first() // .first converts it to a value channel so it can be reused for multiple reports.
    )

    emit:
    multiqc_report = MULTIQC.out.report
}
