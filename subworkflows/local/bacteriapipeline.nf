
include { MAKE_REFERENCE_INDEX } from './make_ref_index'  // this is being called from subworflows




workflow BACTERIAPIPELINE {

    take:
    input // channel: [ tuple val(meta), env(TAXON), path(hits), path(reads), path(reference) ]

    main:
    ch_taxon     = input.map { [it[0], it[1]] } 
    ch_bbsketch  = input.map { [it[0], it[2]] }
    ch_reads     = input.map { [it[0], it[3]] }
    ch_reference = input.map { [it[0], it[4]] }
    ch_versions = Channel.empty()

    MAKE_REFERENCE_INDEX ( ch_reference )
    
    ALIGN_READS_TO_REF (
        ch_reads
        .join(ch_reference)
        .join(MAKE_REFERENCE_INDEX.out.samtools_fai)
        .join(MAKE_REFERENCE_INDEX.out.bwa_index)
    )

    emit:
    picard_dict  = MAKE_REFERENCE_INDEX.out.picard_dict
    samtools_fai = MAKE_REFERENCE_INDEX.out.samtools_fai
    samtools_gzi = MAKE_REFERENCE_INDEX.out.samtools_gzi
    bwa_index    = MAKE_REFERENCE_INDEX.out.bwa_index 
    versions = ch_versions                           // channel: [ versions.yml ]

}

