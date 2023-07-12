include { FASTP           } from '../../modules/nf-core/fastp/main'
include { SPADES          } from '../../modules/nf-core/spades/main'
include { FILTER_ASSEMBLY } from '../../modules/local/filter_assembly'


workflow GENOME_ASSEMBLY {

    take:
    ch_input // channel: [ val(meta), [fastq_1, fastq_2 ]

    main:

    ch_versions = Channel.empty()

    FASTP ( ch_input, [], false, false )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    SPADES (
        FASTP.out.reads.map { [it[0], it[1], [], []] }, // tuple val(meta), path(illumina), path(pacbio), path(nanopore)
        [], // val yml
        []  // val hmm
    )
    ch_versions = ch_versions.mix(SPADES.out.versions.first())

    FILTER_ASSEMBLY (
        SPADES.out.scaffolds
    )
    ch_versions = ch_versions.mix(FILTER_ASSEMBLY.out.versions.first())

    emit:
    reads    = FASTP.out.reads           // channel: [ val(meta), [ fastq_1, fastq_2 ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

