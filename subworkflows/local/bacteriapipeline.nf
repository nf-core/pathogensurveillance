
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

    emit:
    output = MAKE_REFERENCE_INDEX.out.reference_dict


}

