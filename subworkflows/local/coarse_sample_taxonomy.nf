include { BBMAP_SENDSKETCH       } from '../../modules/local/sendsketch'
include { INITIALCLASSIFICATION  } from '../../modules/local/initialclassification'

workflow COARSE_SAMPLE_TAXONOMY {

    take:
    ch_reads_ref  // channel: [ val(meta), [ file(reads) ], val(ref_meta), file(reference) ]

    main:
    ch_versions = Channel.empty()

    ch_fastq = ch_reads_ref.map { it[0..1] }
    BBMAP_SENDSKETCH ( ch_fastq )
    ch_versions = ch_versions.mix(BBMAP_SENDSKETCH.out.versions.toSortedList().map{it[0]})

    INITIALCLASSIFICATION ( BBMAP_SENDSKETCH.out.hits )
    ch_versions = ch_versions.mix(INITIALCLASSIFICATION.out.versions.toSortedList().map{it[0]})

    emit:
    taxon           = INITIALCLASSIFICATION.out.taxon           // channel: [ val(meta), val(taxon) ]
    classification  = INITIALCLASSIFICATION.out.classification  // channel: [ val(meta), val(classification) ]
    hits            = BBMAP_SENDSKETCH.out.hits                 // channel: [ val(meta), file(hits) ]
    versions        = ch_versions                               // channel: [ versions.yml ]
}
