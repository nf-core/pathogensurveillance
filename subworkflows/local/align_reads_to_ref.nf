include { BWA_MEM                            } from '../../modules/nf-core/bwa/mem/main'
include { PICARD_ADDORREPLACEREADGROUPS      } from '../../modules/nf-core/picard/addorreplacereadgroups/main'
include { PICARD_SORTSAM as PICARD_SORTSAM_1 } from '../../modules/nf-core/picard/sortsam/main'
include { PICARD_SORTSAM as PICARD_SORTSAM_2 } from '../../modules/nf-core/picard/sortsam/main'
include { PICARD_MARKDUPLICATES              } from '../../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_INDEX                     } from '../../modules/nf-core/samtools/index/main'

workflow ALIGN_READS_TO_REF {

    take:
    ch_input     // channel: [ val(meta), [ fastq_1, fastq_2 ], ref_meta, reference, reference_index, bam_index ]

    main:

    ch_versions = Channel.empty()

    ch_reads     = ch_input.map { [it[0], it[1]] }
    ch_bwa_index = ch_input.map { [it[0], it[5]] }
    BWA_MEM ( ch_reads, ch_bwa_index, false )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    PICARD_ADDORREPLACEREADGROUPS ( BWA_MEM.out.bam )
    ch_versions = ch_versions.mix(PICARD_ADDORREPLACEREADGROUPS.out.versions.first())

    PICARD_SORTSAM_1 ( PICARD_ADDORREPLACEREADGROUPS.out.bam, 'coordinate' )
    ch_versions = ch_versions.mix(PICARD_SORTSAM_1.out.versions.first())
    
    ch_reference = ch_input.map { it[3] }
    ch_ref_index = ch_input.map { it[4] }
    PICARD_MARKDUPLICATES (
        PICARD_SORTSAM_1.out.bam,
        ch_reference,
        ch_ref_index
    )
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions.first())

    PICARD_SORTSAM_2 ( PICARD_MARKDUPLICATES.out.bam, 'coordinate' )

    SAMTOOLS_INDEX ( PICARD_SORTSAM_2.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    bam      = PICARD_SORTSAM_2.out.bam        // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai   // channel: [ val(meta), [ bai ] ]
    versions = ch_versions                     // channel: [ versions.yml ]
}
