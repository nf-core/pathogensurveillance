include { BBMAP_SENDSKETCH       } from '../../modules/local/sendsketch'
include { INITIALCLASSIFICATION  } from '../../modules/local/initialclassification'

workflow COARSE_SAMPLE_TAXONOMY {

    take:
    ch_reads  // channel: [ val(meta), [ file(reads) ] ]

    main:
    ch_versions = Channel.empty()

    BBMAP_SENDSKETCH ( ch_reads )
    ch_versions = ch_versions.mix(BBMAP_SENDSKETCH.out.versions.toSortedList().map{it[0]})
    ch_depth = BBMAP_SENDSKETCH.out.hits
        .map { [it[0], it[2]] }
    ch_hits = BBMAP_SENDSKETCH.out.hits                                        
        .map { [it[0], it[1]] }                                                 

    INITIALCLASSIFICATION ( ch_hits )
    ch_versions = ch_versions.mix(INITIALCLASSIFICATION.out.versions.toSortedList().map{it[0]})

    emit:
    species         = INITIALCLASSIFICATION.out.species         // channel: [ val(meta), file(taxon) ]
    genera          = INITIALCLASSIFICATION.out.genera          // channel: [ val(meta), file(taxon) ]
    families        = INITIALCLASSIFICATION.out.families        // channel: [ val(meta), file(taxon) ]
    classification  = INITIALCLASSIFICATION.out.classification  // channel: [ val(meta), val(classification) ]
    kingdom         = INITIALCLASSIFICATION.out.kingdom         // channel: [ val(meta), val(kingdom) ]
    hits            = ch_hits                                   // channel: [ val(meta), file(hits) ]
    depth           = ch_depth                                  // channel: [ val(meta), file(depth) ]
    versions        = ch_versions                               // channel: [ versions.yml ]
}
