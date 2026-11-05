/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

// Subworkflow_MakeReferenceIndex
include { PICARD_CREATESEQUENCEDICTIONARY } from '../../../modules/nf-core/picard/createsequencedictionary'
include { SAMTOOLS_FAIDX                  } from '../../../modules/nf-core/samtools/faidx'
include { BWA_INDEX                       } from '../../../modules/nf-core/bwa/index'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Index the reference genome REF.fna
//

workflow REFERENCE_INDEX {
    take:
    reference     // [ val(ref_meta), file(reference) ]

    main:
    ch_versions = channel.empty()

    PICARD_CREATESEQUENCEDICTIONARY ( reference )
    ch_versions = ch_versions.mix (PICARD_CREATESEQUENCEDICTIONARY.out.versions)

    SAMTOOLS_FAIDX ( reference, [[], []], false )

    BWA_INDEX ( reference )

    emit:
    picard_dict   = PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict
    samtools_fai  = SAMTOOLS_FAIDX.out.fai
    samtools_gzi  = SAMTOOLS_FAIDX.out.gzi
    bwa_index     = BWA_INDEX.out.index
    versions      = ch_versions                        // channel: [ versions.yml ]

}
