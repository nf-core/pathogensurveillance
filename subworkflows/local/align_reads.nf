include { BWA_MEM         } from '../../modules/nf-core/bwa/mem/main'
include { PICARD_FORMAT   } from '../../modules/local/picard_format.nf'
include { SAMTOOLS_INDEX  } from '../../modules/nf-core/samtools/index/main'

workflow ALIGN_READS {

    take:
    ch_input     // channel: [ val(meta), [ fastq_1, fastq_2 ], ref_meta, reference, reference_index, bam_index ]

    main:

    versions = Channel.empty()

    samp_ref_combo = ch_input
        .map { [[id: "${it[2].id}_${it[0].id}", ref: it[2], sample: it[0]]] + it } // make composite ID for read/ref combos
    ch_reads = samp_ref_combo
        .map { [it[0], it[2] instanceof Collection && it[2].size() > 2 ? it[2][0..1] : it[2]] } // NOTE: not sure why there are sometimes more than 2 read files. Might need to change once long reads are fullly supported
    ch_bwa_index = samp_ref_combo.map { [it[0], it[6]] }
    BWA_MEM ( ch_reads, ch_bwa_index, false )
    versions = versions.mix(BWA_MEM.out.versions)

    ch_reference = samp_ref_combo.map { [it[0], it[4]] } // channel: [ val(ref_samp_meta), file(reference) ]
    ch_ref_index = samp_ref_combo.map { [it[0], it[5]] } // channel: [ val(ref_samp_meta), file(ref_index) ]
    picard_input = BWA_MEM.out.bam // joined to associated right reference with each sample
        .join(ch_reference)
        .join(ch_ref_index)
    PICARD_FORMAT ( picard_input )

    SAMTOOLS_INDEX ( PICARD_FORMAT.out.bam )
    versions = versions.mix(SAMTOOLS_INDEX.out.versions)

    // Revet combined metas back to seperate ones for sample and reference
    out_bam = PICARD_FORMAT.out.bam        // channel: [ val(ref_samp_meta), [ bam ] ]
        .map { [it[0].sample, it[0].ref, it[1]] }
    out_bai = SAMTOOLS_INDEX.out.bai        // channel: [ val(ref_samp_meta), [ bai ] ]
        .map { [it[0].sample, it[0].ref, it[1]] }

    emit:
    bam      = out_bam        // channel: [ val(meta), val(ref_meta), [ bam ] ]
    bai      = out_bai        // channel: [ val(meta), val(ref_meta), [ bai ] ]
    versions = versions    // channel: [ versions.yml ]
}

