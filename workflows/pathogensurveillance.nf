/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowPathogensurveillance.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified.' }
if (!params.bakta_db && !params.download_bakta_db ) {
    exit 1, "No bakta database specified. Use either '--bakta_db' to point to a local bakta database or use '--download_bakta_db true' to download the Bakta database."
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK              } from '../subworkflows/local/input_check'
include { COARSE_SAMPLE_TAXONOMY   } from '../subworkflows/local/coarse_sample_taxonomy'
include { CORE_GENOME_PHYLOGENY    } from '../subworkflows/local/core_genome_phylogeny'
include { VARIANT_ANALYSIS         } from '../subworkflows/local/variant_analysis'
include { DOWNLOAD_REFERENCES      } from '../subworkflows/local/download_references'
include { ASSIGN_REFERENCES        } from '../subworkflows/local/assign_references'
include { GENOME_ASSEMBLY          } from '../subworkflows/local/genome_assembly'
include { BUSCO_PHYLOGENY          } from '../subworkflows/local/busco_phylogeny'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

include { SRATOOLS_FASTERQDUMP        } from '../modules/local/fasterqdump'
include { MAIN_REPORT                 } from '../modules/local/main_report'
include { RECORD_MESSAGES             } from '../modules/local/record_messages'
include { DOWNLOAD_ASSEMBLIES         } from '../modules/local/download_assemblies'
include { PREPARE_REPORT_INPUT        } from '../modules/local/prepare_report_input'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow PATHOGENSURVEILLANCE {

    ch_versions = Channel.empty()
    messages = Channel.empty()

    // Read in samplesheet, validate and stage input files
    INPUT_CHECK (
        ch_input
    )

    // Download FASTQ files if SRA accessions is provided
    ch_sra = INPUT_CHECK.out.sample_data // meta, [shortread], nanopore, pacbio, sra, ref_meta, reference, reference_refseq, group
        .map { [it[0], it[4]] } // meta, sra
        .filter { it[1] != null }
        .unique()
    SRATOOLS_FASTERQDUMP ( ch_sra )

    // Download reference files if an accession is provided instead of a file
    ch_ref_accessions = INPUT_CHECK.out.sample_data // meta, [shortread], nanopore, pacbio, sra, ref_meta, reference, reference_refseq, group
        .filter { it[7] != null }
        .map { [it[5], it[7]] } // ref_meta, reference_refseq
        .unique()
    DOWNLOAD_ASSEMBLIES ( ch_ref_accessions )
    ch_downloaded_refs = DOWNLOAD_ASSEMBLIES.out.sequence // val(ref_meta), file(downloaded_ref)
        .combine( INPUT_CHECK.out.sample_data.map { [it[5], it[0]] }, by: 0) // val(ref_meta), file(downloaded_ref), val(meta)
        .map { [it[2], it[1]] } // val(meta), file(downloaded_ref)

    // Replace NCBI SRAs/Assembly accessions with downloaded reads and references
    ch_input_parsed = INPUT_CHECK.out.sample_data // meta, [shortread], nanopore, pacbio, sra, ref_meta, reference, reference_refseq, group
        .join(SRATOOLS_FASTERQDUMP.out.reads, remainder:true) // meta, [shortread], nanopore, pacbio, sra, ref_meta, reference, reference_refseq, group, [sra_fastq]
        .join(ch_downloaded_refs, remainder:true) // meta, [shortread], nanopore, pacbio, sra, ref_meta, reference, reference_refseq, group, [sra_fastq], downloaded_ref
        .map { [it[0], it[9] ?: it[1], it[5], it[10] ?: it[6], it[8]] } // meta, [fastq], ref_meta, reference, group_meta

    ch_reads = ch_input_parsed // [val(meta), [file(fastq)], val(ref_meta), file(reference), val(group_meta)]
        .map { it[0..1] }
        .unique()
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // Run FastQC
    FASTQC (
        ch_reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.toSortedList().map{it[0]})

    // Make initial taxonomic classification to decide how to treat sample
    COARSE_SAMPLE_TAXONOMY (
        ch_reads
    )
    ch_versions = ch_versions.mix(COARSE_SAMPLE_TAXONOMY.out.versions)

    // Search for and download reference assemblies for all samples
    DOWNLOAD_REFERENCES (
        COARSE_SAMPLE_TAXONOMY.out.species,
        COARSE_SAMPLE_TAXONOMY.out.genera,
        COARSE_SAMPLE_TAXONOMY.out.families,
        INPUT_CHECK.out.sample_data
    )
    ch_versions = ch_versions.mix(DOWNLOAD_REFERENCES.out.versions)

    // Assign closest reference for samples without a user-assigned reference
    ASSIGN_REFERENCES (
        ch_input_parsed,
        DOWNLOAD_REFERENCES.out.assem_samp_combos,
        DOWNLOAD_REFERENCES.out.sequence,
        DOWNLOAD_REFERENCES.out.signatures,
        COARSE_SAMPLE_TAXONOMY.out.depth
    )
    ch_versions = ch_versions.mix(ASSIGN_REFERENCES.out.versions)
    messages = messages.mix(ASSIGN_REFERENCES.out.messages)

    // Call variants and create SNP-tree and minimum spanning nextwork
    VARIANT_ANALYSIS (
        ASSIGN_REFERENCES.out.sample_data,
        INPUT_CHECK.out.csv
    )
    ch_versions = ch_versions.mix(VARIANT_ANALYSIS.out.versions)
    messages = messages.mix(VARIANT_ANALYSIS.out.messages)

    // Assemble and annotate bacterial genomes
    GENOME_ASSEMBLY (
        ASSIGN_REFERENCES.out.sample_data
            .combine(COARSE_SAMPLE_TAXONOMY.out.kingdom, by: 0)
            .combine(COARSE_SAMPLE_TAXONOMY.out.depth, by: 0)
    )
    ch_versions = ch_versions.mix(GENOME_ASSEMBLY.out.versions)

    // Create core gene phylogeny for bacterial samples
    ref_gffs = DOWNLOAD_REFERENCES.out.assem_samp_combos
        .combine(DOWNLOAD_REFERENCES.out.gff, by: 0) // [ val(ref_meta), val(meta), file(gff) ]
        .map { it[1..2] } // [ val(meta), file(gff) ]
        .groupTuple() // [ val(meta), [file(gff)] ]
    gff_and_group = ASSIGN_REFERENCES.out.sample_data  // [val(meta), [file(fastq)], val(ref_meta), file(reference), val(group_meta)]
        .combine(GENOME_ASSEMBLY.out.gff, by: 0) // [val(meta), [file(fastq)], val(ref_meta), file(reference), val(group_meta), file(gff)]
        .combine(ref_gffs, by: 0) // [val(meta), [file(fastq)], val(ref_meta), file(reference), val(group_meta), file(gff), [file(ref_gff)] ]
        .combine(COARSE_SAMPLE_TAXONOMY.out.depth, by:0) // [val(meta), [file(fastq)], val(ref_meta), file(reference), val(group_meta), file(gff), [file(ref_gff)], val(depth)]
        .map { [it[0], it[5], it[4], it[6], it[7]] } // [ val(meta), file(gff), val(group_meta), [file(ref_gff)], val(depth) ]
    CORE_GENOME_PHYLOGENY (
        gff_and_group,
        INPUT_CHECK.out.csv
    )
    ch_versions = ch_versions.mix(CORE_GENOME_PHYLOGENY.out.versions)
    messages  = messages.mix(CORE_GENOME_PHYLOGENY.out.messages)

    // Read2tree BUSCO phylogeny for eukaryotes
    ref_metas = DOWNLOAD_REFERENCES.out.assem_samp_combos // val(ref_meta), val(meta)
        .map { [it[1], it[0]] } // val(meta), val(ref_meta)
        .groupTuple() // val(meta), [val(ref_meta)]
    busco_input = ASSIGN_REFERENCES.out.sample_data  // val(meta), [file(fastq)], val(ref_meta), file(reference), val(group_meta)
        .combine(ref_metas, by: 0) // val(meta), [file(fastq)], val(ref_meta), file(reference), val(group_meta), [val(ref_meta)]
        .combine(COARSE_SAMPLE_TAXONOMY.out.kingdom, by: 0) // val(meta), [file(fastq)], val(ref_meta), file(reference), val(group_meta), [val(ref_meta)], val(kingdom)
        .combine(COARSE_SAMPLE_TAXONOMY.out.depth, by:0) // val(meta), [file(fastq)], val(ref_meta), file(reference), val(group_meta), [val(ref_meta)], val(kingdom), val(depth)
        .map { it[0..1] + it[4..7] } // val(meta), [file(fastq)], val(group_meta), [val(ref_meta)], val(kingdom), val(depth)
    //BUSCO_PHYLOGENY (
    //    busco_input,
    //    DOWNLOAD_REFERENCES.out.sequence
    //)

    // Save version info
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )


    // MultiQC
    workflow_summary    = WorkflowPathogensurveillance.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowPathogensurveillance.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.collect().ifEmpty([]),
        ch_multiqc_custom_config.collect().ifEmpty([]),
        ch_multiqc_logo.collect().ifEmpty([])
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)

    // Save error/waring/message info
    RECORD_MESSAGES (
        messages.collect(sort:true, flat:false)
    )

    // Create main summary report
    report_samp_data = ASSIGN_REFERENCES.out.sample_data // meta, fastq, ref_meta, reference, group_meta
        .combine(COARSE_SAMPLE_TAXONOMY.out.hits, by:0) // meta, fastq, ref_meta, reference, group_meta, sendsketch
        .combine(DOWNLOAD_REFERENCES.out.stats, by:0) // meta, fastq, ref_meta, reference, group_meta, sendsketch, ref_stats
        .map { [it[2]] + it[0..1] + it[3..6] } // ref_meta, meta, fastq, reference, group_meta, sendsketch, ref_stats
        .join(GENOME_ASSEMBLY.out.quast, remainder: true, by:0) // ref_meta, meta, fastq, reference, group_meta, sendsketch, ref_stats, quast
        .map { [it[4]] + it[0..3] + it[5..7]} // group_meta, ref_meta, meta, fastq, reference, sendsketch, ref_stats, quast
        .groupTuple() // group_meta, [ref_meta], [meta], [fastq], [reference], [sendsketch], [ref_stats], [quast]
    report_variant_data = VARIANT_ANALYSIS.out.results // group_meta, ref_meta, vcf, align, tree
        .groupTuple() // group_meta, [ref_meta], [vcf], [align], [tree]
    report_group_data = ASSIGN_REFERENCES.out.ani_matrix // group_meta, ani_matrix
        .join(CORE_GENOME_PHYLOGENY.out.phylogeny, remainder:true) // group_meta, ani_matrix, core_phylo
        .join(CORE_GENOME_PHYLOGENY.out.pocp, remainder:true) // group_meta, ani_matrix, core_phylo, popc
        .join(ASSIGN_REFERENCES.out.assigned_refs, remainder:true) // group_meta, ani_matrix, core_phylo, pocp, assigned_refs
    report_in = report_samp_data // group_meta, [ref_meta], [meta], [fastq], [reference], [sendsketch], [ref_stats], [quast]
        .join(report_variant_data, remainder: true) // group_meta, [ref_meta], [meta], [fastq], [reference], [sendsketch], [ref_stats], [quast], [ref_meta], [vcf], [align], [tree]
        .map { it.size() == 12 ? it : it[0..7] + [[], [], [], []] } // group_meta, [ref_meta], [meta], [fastq], [reference], [sendsketch], [ref_stats], [quast], [ref_meta], [vcf], [align], [tree]
        .join(report_group_data, remainder: true) // group_meta, [ref_meta], [meta], [fastq], [reference], [sendsketch], [ref_stats], [quast], [ref_meta], [vcf], [align], [tree], ani_matrix, core_phylo, pocp, passigned_refs
        .map { it.size() == 16 ? it : it[0..11] + [[], [], [], []] } // group_meta, [ref_meta], [meta], [fastq], [reference], [sendsketch], [ref_stats], [quast], [ref_meta], [vcf], [align], [tree], ani_matrix, core_phylo, pocp, assigned_refs
        .map { it[0..1] + it[5..7] + it[9..15] } // group_meta, [ref_meta], [sendsketch], [ref_stats], [quast], [vcf], [align], [tree], ani_matrix, core_phylo, pocp, assigned_refs
        .map { [
            it[0],
            it[1].findAll{ it.id != null }.unique(),
            it[2],
            it[3],
            it[4].findAll{ it != null },
            it[5].findAll{ it != null },
            it[6].findAll{ it != null },
            it[7].findAll{ it != null },
            it[8],
            it[9] == null ? [] : it[9],
            it[10],
            it[11]
         ] } // group_meta, [ref_meta], [sendsketch], [ref_stats], [quast], [vcf], [align], [tree], ani_matrix, core_phylo, pocp, assigned_refs

    PREPARE_REPORT_INPUT (
        report_in,
        INPUT_CHECK.out.csv,
        MULTIQC.out.data,
        MULTIQC.out.plots,
        MULTIQC.out.report,
        CUSTOM_DUMPSOFTWAREVERSIONS.out.yml.first(), // .first converts it to a value channel so it can be reused for multiple reports.
        RECORD_MESSAGES.out.tsv
    )

    MAIN_REPORT (
        PREPARE_REPORT_INPUT.out.report_input,
        Channel.fromPath("${projectDir}/assets/main_report", checkIfExists: true).first() // .first converts it to a value channel so it can be reused for multiple reports.
    )

}










/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.adaptivecard(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
