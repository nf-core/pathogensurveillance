include { BBMAP_SENDSKETCH      } from '../../modules/local/sendsketch'
include { INITIALCLASSIFICATION  } from '../../modules/local/initialclassification'

workflow COARSE_SAMPLE_TAXONOMY {

    take:
    ch_fastq  // channel: [ val(meta), [ fastq ] ]

    main:
    ch_versions = Channel.empty()

    BBMAP_SENDSKETCH ( ch_fastq )
    ch_versions = ch_versions.mix(BBMAP_SENDSKETCH.out.versions.first())

    INITIALCLASSIFICATION ( BBMAP_SENDSKETCH.out.hits )
    ch_versions = ch_versions.mix(INITIALCLASSIFICATION.out.versions.first())

    emit:
    taxon             = INITIALCLASSIFICATION.out.taxon  // channel: [ val(meta), [ taxon ] ]
    classification    = INITIALCLASSIFICATION.out.classification  // channel: [ val(meta), [ classification ] ]
    hits              = BBMAP_SENDSKETCH.out.hits     // channel: [ val(meta), [ hits ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

