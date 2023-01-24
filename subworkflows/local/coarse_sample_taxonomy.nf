include { BBMAP_SENDSKETCH      } from '../../modules/local/sendsketch'
// include { INITIAL_CLASSIFICATION  } from '../../modules/local/initial_classification/main'

workflow COARSE_SAMPLE_TAXONOMY {

    take:
    ch_fastq  // channel: [ val(meta), [ fastq ] ]

    main:
    ch_versions = Channel.empty()

    BBMAP_SENDSKETCH ( ch_fastq )
    ch_versions = ch_versions.mix(BBMAP_SENDSKETCH.out.versions.first())

    // INITIAL_CLASSIFICATION ( BBMAP_SENDSKETCH.out.hits )
    // ch_versions = ch_versions.mix(INITIAL_CLASSIFICATION.out.versions.first())

    emit:
    // taxon    = INITIAL_CLASSIFICATION.out.bam  // channel: [ val(meta), [ taxon ] ]
    // class    = INITIAL_CLASSIFICATION.out.bai  // channel: [ val(meta), [ class ] ]
    hits     = BBMAP_SENDSKETCH.out.hits     // channel: [ val(meta), [ hits ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

