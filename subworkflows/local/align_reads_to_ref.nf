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
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.toSortedList().map{it[0]})

    PICARD_ADDORREPLACEREADGROUPS ( BWA_MEM.out.bam )
    ch_versions = ch_versions.mix(PICARD_ADDORREPLACEREADGROUPS.out.versions.toSortedList().map{it[0]})

    PICARD_SORTSAM_1 ( PICARD_ADDORREPLACEREADGROUPS.out.bam, 'coordinate' )
    ch_versions = ch_versions.mix(PICARD_SORTSAM_1.out.versions.toSortedList().map{it[0]})
    
    ch_reference = ch_input.map { [it[0], it[3]] } // channel: [ val(meta), file(reference) ]
    ch_ref_index = ch_input.map { [it[0], it[4]] } // channel: [ val(meta), file(ref_index) ]
    picard_input = PICARD_SORTSAM_1.out.bam // joined to associated right reference with each sample
        .join(ch_reference)
        .join(ch_ref_index)

    PICARD_MARKDUPLICATES (
        picard_input.map { it[0..1] },
        picard_input.map { it[2] },
        picard_input.map { it[3] }
    )
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions.toSortedList().map{it[0]})

    PICARD_SORTSAM_2 ( PICARD_MARKDUPLICATES.out.bam, 'coordinate' )

    SAMTOOLS_INDEX ( PICARD_SORTSAM_2.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.toSortedList().map{it[0]})

    emit:
    bam      = PICARD_SORTSAM_2.out.bam        // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai   // channel: [ val(meta), [ bai ] ]
    versions = ch_versions                     // channel: [ versions.yml ]
}

