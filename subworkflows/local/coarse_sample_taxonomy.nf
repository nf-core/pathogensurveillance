include { BBMAP_SENDSKETCH      } from '../../modules/local/sendsketch'
include { INITIALCLASSIFICATION  } from '../../modules/local/initialclassification'

workflow COARSE_SAMPLE_TAXONOMY {

    take:
    ch_fastq  // channel: [ val(meta), [ fastq ] ]

    main:
    ch_versions = Channel.empty()

    BBMAP_SENDSKETCH ( ch_fastq )
    ch_versions = ch_versions.mix(BBMAP_SENDSKETCH.out.versions.first())
    ch_hits_and_reads = BBMAP_SENDSKETCH.out.hits.join(ch_fastq)

    INITIALCLASSIFICATION ( ch_hits_and_reads )
    ch_versions = ch_versions.mix(INITIALCLASSIFICATION.out.versions.first())

    emit:
    result          = INITIALCLASSIFICATION.out.result  // channel: [ val(meta), env(TAXON), path(hits), path(reads) ]
    classification  = INITIALCLASSIFICATION.out.classification  // channel: [ val(meta), [ classification ] ]
    hits            = BBMAP_SENDSKETCH.out.hits     // channel: [ val(meta), [ hits ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

