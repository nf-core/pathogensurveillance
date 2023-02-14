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
    ch_reference     // [ val(meta), val(ref_meta), file(reference) ]

    main:
    ch_versions = Channel.empty()

    // make channel for each unique reference
    ch_unique_ref = ch_reference
        .map { it[1..2] }
        .groupTuple() // make unique based on reference metadata
        .map { [it[0],it[1].sort()[0]] } // picks first sorted ref link so the cache works
    PICARD_CREATESEQUENCEDICTIONARY ( ch_unique_ref )
    ch_versions = ch_versions.mix (PICARD_CREATESEQUENCEDICTIONARY.out.versions)

    SAMTOOLS_FAIDX ( ch_unique_ref )
    ch_versions = ch_versions.mix (SAMTOOLS_FAIDX.out.versions)

    BWA_INDEX ( ch_unique_ref ) 
    ch_versions = ch_versions.mix (BWA_INDEX.out.versions)

    // duplicate and recombine unique reference results with sample data
    //    odd cross + map needed since join does not duplicate things
    picard_dict = PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict
        .cross(ch_reference.map { [it[1], it[0]] } )
        .map { [it[1][1], it[0][1]] }
    samtools_fai = SAMTOOLS_FAIDX.out.fai
        .cross(ch_reference.map { [it[1], it[0]] } )
        .map { [it[1][1], it[0][1]] } 
    samtools_gzi = SAMTOOLS_FAIDX.out.gzi
        .cross(ch_reference.map { [it[1], it[0]] } )
        .map { [it[1][1], it[0][1]] } 
    bwa_index = BWA_INDEX.out.index
        .cross(ch_reference.map { [it[1], it[0]] } )
        .map { [it[1][1], it[0][1]] } 

    emit:
    picard_dict
    samtools_fai
    samtools_gzi
    bwa_index
    versions = ch_versions                        // channel: [ versions.yml ]

}
