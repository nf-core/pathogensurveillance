include { KHMER_TRIMLOWABUND       } from '../../modules/local/khmer_trimlowabund'
include { SOURMASH_SKETCH          } from '../../modules/nf-core/sourmash/sketch/main'
include { SOURMASH_COMPARE         } from '../../modules/nf-core/sourmash/compare/main'

workflow ASSIGN_REFERENCES {

    take:
    sample_data  // [val(meta), [file(fastq)], val(ref_meta), file(reference), val(groups)]
    assem_samp_combos // [val(assem_id), val(meta)] for each unique combination
    sequence // [val(assem_meta), file(fna)] for each assembly
    signatures // [val(assem_meta), file(sig)] for each assembly

    main:
    ch_versions = Channel.empty()

    KHMER_TRIMLOWABUND (
        sample_data
            .map { it[0..1] } // [val(meta), [file(fastq)]], possibly duplicated
            .distinct() // [val(meta), [file(fastq)]], one per sample
    )

    SOURMASH_SKETCH (
        KHMER_TRIMLOWABUND.out.sequence
    )

    sample_sigs = sample_data
        .combine(SOURMASH_SKETCH.out.signatures, by: 0)
        .map { [it[4], it[5]] }
        .groupTuple() // group, [sig]

    assem_sigs = assem_samp_combos
        .combine(signatures.map { [it[0].id, it[1]] }, by: 0)
        .map { [it[1], it[0], it[2]] } // meta, assem, sig
        .combine(sample_data, by: 0) // meta, assem, sig, fastq, ref_meta, ref, groups
        .map { [it[6], it[2]] } // group, sig
        .groupTuple() // group, [sig]
        .map { [it[0], it[1].unique()] }

    group_sigs = sample_sigs
        .join(assem_sigs)
        .map { [[id: it[0]], it[1] + it[2]] }

    SOURMASH_COMPARE (
        group_sigs,
        [], // file_list (optional)
        true, // save numpy matrix
        true  // save CSV
    )


    emit:
    versions        = ch_versions                      // channel: [ versions.yml ]
}
