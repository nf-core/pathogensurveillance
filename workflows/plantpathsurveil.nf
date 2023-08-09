/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowPlantpathsurveil.initialise(params, log)

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
include { VARIANT_CALLING_ANALYSIS } from '../subworkflows/local/variant_calling_analysis'
include { DOWNLOAD_REFERENCES      } from '../subworkflows/local/download_references'
include { ASSIGN_REFERENCES        } from '../subworkflows/local/assign_references'
include { GENOME_ASSEMBLY          } from '../subworkflows/local/genome_assembly'

include { MAIN_REPORT              } from '../modules/local/mainreport'

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

workflow PLANTPATHSURVEIL {

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

    // Assign closest reference for samples without a user-assigned reference
    ASSIGN_REFERENCES (
        INPUT_CHECK.out.sample_data,
        DOWNLOAD_REFERENCES.out.assem_samp_combos,
        DOWNLOAD_REFERENCES.out.sequence,
        DOWNLOAD_REFERENCES.out.signatures
    )

    // Call variants and create SNP-tree and minimum spanning nextwork
    VARIANT_CALLING_ANALYSIS (
        ASSIGN_REFERENCES.out.sample_data,
        ch_input
    )

    // Assemble and annotate bacterial genomes
    GENOME_ASSEMBLY (                                                           
        ASSIGN_REFERENCES.out.sample_data                           
            .combine(COARSE_SAMPLE_TAXONOMY.out.kingdom, by: 0)                                                   
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
        .map { [it[0], it[5], it[4], it[6]] } // [ val(meta), file(gff), val(group_meta), [file(ref_gff)] ]            
    CORE_GENOME_PHYLOGENY (                                                     
        gff_and_group,                             
        ch_input                                                          
    )                                                                           

    // Read2tree phylogeny for eukaryotes
    //READ2TREE_ANALYSIS (
    //)

    // Create main summary report
    report_in = VARIANT_CALLING_ANALYSIS.out.phylogeny // [ group_meta, ref_meta, tree ]
        .groupTuple() // [ group_meta, [ref_meta], [tree] ]
        .join(ASSIGN_REFERENCES.out.ani_matrix) // [ group_meta, [ref_meta], [tree], ani_matrix ]
        .join(CORE_GENOME_PHYLOGENY.out.phylogeny, remainder: true) // [ group_meta, [ref_meta], [snp_tree], ani_matrix, core_tree ]

    MAIN_REPORT ( 
        report_in,
        ch_input,
        DOWNLOAD_REFERENCES.out.stats
    )  
    
    // Save version info
    CUSTOM_DUMPSOFTWAREVERSIONS (                                               
        ch_versions.unique().collect(sort:true)
    )
                                                                          
    // MultiQC
    //workflow_summary    = WorkflowPlantpathsurveil.paramsSummaryMultiqc(workflow, summary_params)
    //ch_workflow_summary = Channel.value(workflow_summary)

    //methods_description    = WorkflowPlantpathsurveil.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    //ch_methods_description = Channel.value(methods_description)

    //ch_multiqc_files = Channel.empty()
    //ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    //ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    //ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    //MULTIQC (
    //    ch_multiqc_files.collect(),
    //    ch_multiqc_config.collect().ifEmpty([]),
    //    ch_multiqc_custom_config.collect().ifEmpty([]),
    //    ch_multiqc_logo.collect().ifEmpty([])
    //)
    //multiqc_report = MULTIQC.out.report.toList()
    //ch_versions    = ch_versions.mix(MULTIQC.out.versions)


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
