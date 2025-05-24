include { BWA_MEM         } from '../../../modules/nf-core/bwa/mem'
include { PICARD_FORMAT   } from '../../../modules/local/picard_format'
include { SAMTOOLS_INDEX  } from '../../../modules/nf-core/samtools/index'

workflow ALIGN_READS {

    take:
    ch_input // channel: [ val(meta), [ fastq_1, fastq_2 ], ref_meta, reference, reference_index, bam_index ]

    main:

    versions = Channel.empty()

    // Addd composite ID for read/ref combos to input
    samp_ref_combo = ch_input
        .map { meta, reads, ref_meta, ref, ref_index, bam_index ->
            [[id: "${ref_meta.id}--${meta.id}", ref: ref_meta, sample: meta], meta, reads, ref_meta, ref, ref_index, bam_index]
        }

    // Make channel with reads, removing any extra fastq files, which are ususally unpaired reads
    ch_reads = samp_ref_combo
        .map { combined_meta, meta, fastqs, ref_meta, reference, ref_index, bam_index ->
            [combined_meta, fastqs instanceof Collection && fastqs.size() > 2 ? fastqs[0..1] : fastqs]
        }

    // Align reads with BWA
    ch_bwa_index = samp_ref_combo.map { combined_meta, meta, fastqs, ref_meta, reference, ref_index, bam_index ->
        [combined_meta, bam_index]
    }
    BWA_MEM ( ch_reads, ch_bwa_index, [[], []], false )
    versions = versions.mix(BWA_MEM.out.versions)

    // Run a series of picard commands to sort and filter variants
    ch_reference = samp_ref_combo.map { combined_meta, meta, fastqs, ref_meta, reference, ref_index, bam_index ->
        [combined_meta, reference]
    }
    ch_ref_index = samp_ref_combo.map { combined_meta, meta, fastqs, ref_meta, reference, ref_index, bam_index ->
        [combined_meta, ref_index]
    }
    picard_input = BWA_MEM.out.bam
        .join(ch_reference)
        .join(ch_ref_index)
    PICARD_FORMAT ( picard_input )
    versions = versions.mix(PICARD_FORMAT.out.versions)

    SAMTOOLS_INDEX ( PICARD_FORMAT.out.bam )
    versions = versions.mix(SAMTOOLS_INDEX.out.versions)

    // Revet combined metas back to seperate ones for sample and reference
    out_bam = PICARD_FORMAT.out.bam
        .map { combined_meta, bam ->
            [combined_meta.sample, combined_meta.ref, bam]
        }
    out_csi = SAMTOOLS_INDEX.out.csi
        .map { combined_meta, csi ->
            [combined_meta.sample, combined_meta.ref, csi]
        }

    emit:
    bam      = out_bam        // channel: [ val(meta), val(ref_meta), [ bam ] ]
    csi      = out_csi        // channel: [ val(meta), val(ref_meta), [ csi ] ]
    versions = versions    // channel: [ versions.yml ]
}

