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
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }



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

include { MAIN_REPORT as MAIN_REPORT_1 } from '../modules/local/main_report'
include { MAIN_REPORT as MAIN_REPORT_2 } from '../modules/local/main_report'

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


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow PATHOGENSURVEILLANCE {

    ch_versions = Channel.empty()

    // Read in samplesheet, validate and stage input files
    INPUT_CHECK (
        ch_input
    )
    ch_reads = INPUT_CHECK.out.sample_data // [val(meta), [file(fastq)], val(ref_meta), file(reference), val(group_meta)]
        .map { it[0..1] }
        .distinct()
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
        COARSE_SAMPLE_TAXONOMY.out.families
    )
    ch_versions = ch_versions.mix(DOWNLOAD_REFERENCES.out.versions)

    // Create main summary report                                               
    //MAIN_REPORT_1 (                                                               
    //    INPUT_CHECK.out.sample_data.map {[ it[4], it[2], null ]}.groupTuple().map {it + [null, null]}, 
    //    ch_input,                                                               
    //    DOWNLOAD_REFERENCES.out.stats                                           
    //)                                                                           

    // Assign closest reference for samples without a user-assigned reference
    ASSIGN_REFERENCES (
        INPUT_CHECK.out.sample_data,
        DOWNLOAD_REFERENCES.out.assem_samp_combos,
        DOWNLOAD_REFERENCES.out.sequence,
        DOWNLOAD_REFERENCES.out.signatures,
        COARSE_SAMPLE_TAXONOMY.out.depth
    )
    ch_versions = ch_versions.mix(ASSIGN_REFERENCES.out.versions)

    // Call variants and create SNP-tree and minimum spanning nextwork
    VARIANT_ANALYSIS (
        ASSIGN_REFERENCES.out.sample_data,
        ch_input
    )
    ch_versions = ch_versions.mix(VARIANT_ANALYSIS.out.versions)

    // Assemble and annotate bacterial genomes
    GENOME_ASSEMBLY (                                                           
        ASSIGN_REFERENCES.out.sample_data                           
            .combine(COARSE_SAMPLE_TAXONOMY.out.kingdom, by: 0)
            .combine(COARSE_SAMPLE_TAXONOMY.out.depth, by: 0)                                                   
    )                                                                           
    ch_versions = ch_versions.mix(GENOME_ASSEMBLY.out.versions)                 

    // Create core gene phylogeny for bacterial samples
    ref_gffs = DOWNLOAD_REFERENCES.out.assem_samp_combos
        .combine(DOWNLOAD_REFERENCES.out.gff, by: 0) // [ val(genome_id), val(meta), file(gff) ]
        .map { it[1..2] } // [ val(meta), file(gff) ]
        .groupTuple() // [ val(meta), [file(gff)] ]
    gff_and_group = ASSIGN_REFERENCES.out.sample_data  // [val(meta), [file(fastq)], val(ref_meta), file(reference), val(group_meta)]
        .combine(GENOME_ASSEMBLY.out.gff, by: 0) // [val(meta), [file(fastq)], val(ref_meta), file(reference), val(group_meta), file(gff)]
        .combine(ref_gffs, by: 0) // [val(meta), [file(fastq)], val(ref_meta), file(reference), val(group_meta), file(gff), [file(ref_gff)] ]
        .combine(COARSE_SAMPLE_TAXONOMY.out.depth, by:0) // [val(meta), [file(fastq)], val(ref_meta), file(reference), val(group_meta), file(gff), [file(ref_gff)], val(depth)]
        .map { [it[0], it[5], it[4], it[6], it[7]] } // [ val(meta), file(gff), val(group_meta), [file(ref_gff)], val(depth) ]            
    CORE_GENOME_PHYLOGENY (                                                     
        gff_and_group,                             
        ch_input                                                          
    )
    ch_versions = ch_versions.mix(CORE_GENOME_PHYLOGENY.out.versions)

    // Read2tree phylogeny for eukaryotes
    //READ2TREE_ANALYSIS (
    //)

    // Save version info
    CUSTOM_DUMPSOFTWAREVERSIONS (                                               
        ch_versions.unique().collect(sort:true)
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

    // Create main summary report                                               
    grouped_sendsketch = INPUT_CHECK.out.sample_data // meta, fastq, ref_meta, reference, group_meta
        .combine(COARSE_SAMPLE_TAXONOMY.out.hits, by:0) // meta, fastq, ref_meta, reference, group_meta, sendsketch
        .map { it[4..5] } // group_meta, sendsketch                             
        .groupTuple() // group_meta, [sendsketch]
    report_in = VARIANT_ANALYSIS.out.phylogeny // group_meta, ref_meta, tree    
        .combine(VARIANT_ANALYSIS.out.snp_align, by:0..1) // group_meta, ref_meta, tree, snp_align
        .combine(VARIANT_ANALYSIS.out.vcf, by:0..1) // group_meta, ref_meta, tree, snp_align, vcf
        .map { [it[1], it[0]] + it[2..4] } // ref_meta, group_meta, tree, snp_align, vcf
        .join(GENOME_ASSEMBLY.out.quast, remainder: true) // ref_meta, group_meta, tree, snp_align, vcf, quast
        .map { [it[1], it[0]] + it[2..5]}
        .groupTuple() // group_meta, [ref_meta], [tree], [snp_align], [vcf], [quast]
        .join(ASSIGN_REFERENCES.out.ani_matrix, remainder: true) // group_meta, [ref_meta], [tree], [snp_align], [vcf], [quast] ani_matrix
        .join(CORE_GENOME_PHYLOGENY.out.phylogeny, remainder: true) // group_meta, [ref_meta], [snp_tree], [snp_align], [vcf], [quast], ani_matrix, core_tree
        .join(grouped_sendsketch, remainder: true) // group_meta, [ref_meta], [snp_tree], [snp_align], [vcf], [quast], ani_matrix, core_tree, [sendsketch]
                                                                     
    MAIN_REPORT_2 (                                                             
        report_in,                                                              
        ch_input,                                                               
        DOWNLOAD_REFERENCES.out.stats,                                          
        MULTIQC.out.data,                                                       
        MULTIQC.out.plots,                                                      
        MULTIQC.out.report,
        CUSTOM_DUMPSOFTWAREVERSIONS.out.yml
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
