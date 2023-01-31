/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

// Subworkflow_MakeReferenceIndex
include { PICARD_CREATESEQUENCEDICTIONARY } from '../../modules/nf-core/picard/createsequencedictionary/main'
include { SAMTOOLS_FAIDX                  } from '../../modules/nf-core/samtools/faidx/main'
include { BWA_INDEX                       } from '../../modules/nf-core/bwa/index/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Index the reference genome REF.fna
//

workflow MAKE_REFERENCE_INDEX{
    take:
    ch_reference     // tuple val(meta), path(fasta)

    main:
    ch_versions = Channel.empty()

    //
    // Rename reference to REF.fna and format for BWA
    //


    PICARD_CREATESEQUENCEDICTIONARY ( ch_reference )
    ch_versions = ch_versions.mix (PICARD_CREATESEQUENCEDICTIONARY.out.versions)

    SAMTOOLS_FAIDX ( ch_reference )
    ch_versions = ch_versions.mix (SAMTOOLS_FAIDX.out.versions)

    BWA_INDEX ( ch_reference ) 
    ch_versions = ch_versions.mix (BWA_INDEX.out.versions)

    emit:
    reference_dict        = PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict



    versions = ch_versions                        // channel: [ versions.yml ]

}
