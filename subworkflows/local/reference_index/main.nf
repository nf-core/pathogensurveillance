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
include { BWAMEM3_INDEX                   } from '../../../modules/nf-core/bwamem3/index/main'


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

    PICARD_CREATESEQUENCEDICTIONARY ( reference )

    SAMTOOLS_FAIDX ( reference.map {meta, fasta -> [meta, fasta, []]}, false )

    if (params.aligner == 'bwamem3') {
        BWAMEM3_INDEX ( reference )
        aligner_index = BWAMEM3_INDEX.out.index
    } else {
        BWA_INDEX ( reference )
        aligner_index = BWA_INDEX.out.index
    }

    emit:
    picard_dict   = PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict
    samtools_fai  = SAMTOOLS_FAIDX.out.fai
    samtools_gzi  = SAMTOOLS_FAIDX.out.gzi
    bwa_index     = aligner_index

}
