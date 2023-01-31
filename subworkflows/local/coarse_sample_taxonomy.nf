include { BBMAP_SENDSKETCH      } from '../../modules/local/sendsketch'
include { INITIALCLASSIFICATION  } from '../../modules/local/initialclassification'

workflow COARSE_SAMPLE_TAXONOMY {

    take:
    ch_reads_ref  // channel: [ val(meta), [ fastq ], reference ]

    main:
    ch_versions = Channel.empty()

    ch_fastq = ch_reads_ref.map { it[0..1] }
    BBMAP_SENDSKETCH ( ch_fastq )
    ch_versions = ch_versions.mix(BBMAP_SENDSKETCH.out.versions.first())
    ch_hits_reads_ref = BBMAP_SENDSKETCH.out.hits.join(ch_reads_ref)

    INITIALCLASSIFICATION ( ch_hits_reads_ref )
    ch_versions = ch_versions.mix(INITIALCLASSIFICATION.out.versions.first())

    emit:
    result          = INITIALCLASSIFICATION.out.result  // channel: [ val(meta), env(TAXON), path(hits), path(reads), path(reference) ]
    classification  = INITIALCLASSIFICATION.out.classification  // channel: [ val(meta), [ classification ] ]
    hits            = BBMAP_SENDSKETCH.out.hits     // channel: [ val(meta), [ hits ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

