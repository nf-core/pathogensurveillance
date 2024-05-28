include { BBMAP_SENDSKETCH       } from '../../modules/local/sendsketch'
include { INITIAL_CLASSIFICATION } from '../../modules/local/initial_classification'

workflow COARSE_SAMPLE_TAXONOMY {

    take:
    ch_reads  // meta, [reads]

    main:
    ch_versions = Channel.empty()

    BBMAP_SENDSKETCH ( ch_reads )
    ch_versions = ch_versions.mix(BBMAP_SENDSKETCH.out.versions.toSortedList().map{it[0]})
    ch_depth = BBMAP_SENDSKETCH.out.hits
        .map { [it[0], it[2]] }
    ch_hits = BBMAP_SENDSKETCH.out.hits
        .map { [it[0], it[1]] }

    INITIAL_CLASSIFICATION ( ch_hits )
    ch_versions = ch_versions.mix(INITIAL_CLASSIFICATION.out.versions.toSortedList().map{it[0]})

    emit:
    species         = INITIAL_CLASSIFICATION.out.species        // val(meta), file(taxon)
    genera          = INITIAL_CLASSIFICATION.out.genera         // val(meta), file(taxon)
    families        = INITIAL_CLASSIFICATION.out.families       // val(meta), file(taxon)
    classification  = INITIAL_CLASSIFICATION.out.classification // val(meta), val(classification)
    kingdom         = INITIAL_CLASSIFICATION.out.kingdom        // val(meta), val(kingdom)
    hits            = ch_hits                                   // val(meta), file(hits)
    depth           = ch_depth                                  // val(meta), file(depth)
    versions        = ch_versions                               // versions.yml
}
