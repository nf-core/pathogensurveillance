include { KHMER_TRIMLOWABUND       } from '../../modules/local/khmer_trimlowabund'
include { ASSIGN_GROUP_REFERENCES  } from '../../modules/local/assign_group_references'
include { SOURMASH_SKETCH as SOURMASH_SKETCH_READS    } from '../../modules/nf-core/sourmash/sketch/main'
include { SOURMASH_SKETCH as SOURMASH_SKETCH_GENOME   } from '../../modules/nf-core/sourmash/sketch/main'
include { SOURMASH_COMPARE         } from '../../modules/local/sourmash_compare'

workflow ASSIGN_REFERENCES {

    take:
    sample_data  // [val(meta), [file(fastq)], val(ref_meta), file(reference), val(group_meta)]
    assem_samp_combos // [val(assem_id), val(meta)] for each unique combination
    sequence // [val(assem_meta), file(fna)] for each assembly
    signatures // [val(assem_meta), file(sig)] for each assembly

    main:
    ch_versions = Channel.empty()

    // Trim rare k-mers from raw reads
    KHMER_TRIMLOWABUND (
        sample_data
            .map { it[0..1] } // [val(meta), [file(fastq)]], possibly duplicated
            .unique() // [val(meta), [file(fastq)]], one per sample
    )

    // Create signature for each sample
    SOURMASH_SKETCH_READS (
        KHMER_TRIMLOWABUND.out.sequence
    )

    // Create signature for each user-defined reference genome
    user_refs = sample_data
        .filter { it[3] != null }
        .map { it[2..3] } // val(ref_meta), file(ref)
        .unique()
    SOURMASH_SKETCH_GENOME (                                                           
        user_refs                                         
    )                                                                           
    
    // Make list of user-defined reference signatures for each group                            
    user_sigs = sample_data
        .map { [it[2], it[4]] } // ref_meta, group_meta                                                  
        .combine(SOURMASH_SKETCH_GENOME.out.signatures, by: 0) // ref_meta, group_meta, sig               
        .map { it[1..2] }  // [group_meta, sig], possibly duplicated
        .unique()                                           
        .groupTuple() // group_meta, [sig]                                      

    // Make list of sample signatures for each group
    sample_sigs = sample_data
        .combine(SOURMASH_SKETCH_READS.out.signatures, by: 0)
        .map { [it[4], it[5]] }
        .groupTuple() // group_meta, [sig]

    // Make list of downloaded reference genome signatures for each group
    assem_sigs = assem_samp_combos
        .combine(signatures.map { [it[0].id, it[1]] }, by: 0)
        .map { [it[1], it[0], it[2]] } // meta, assem, sig
        .combine(sample_data, by: 0) // meta, assem, sig, fastq, ref_meta, ref, group_meta
        .map { [it[6], it[2]] } // group_meta, sig
        .groupTuple() // group_meta, [sig]
        .map { [it[0], it[1].unique()] }

    // Combine all signatures for each group
    group_sigs = sample_sigs
        .join(assem_sigs)
        .join(user_sigs)
        .map { [it[0], it[1] + it[2] + it[3]] }

    // Compare all genomes/samples to eachother to create an ANI matrix
    SOURMASH_COMPARE (
        group_sigs,
        [], // file_list (optional)
        true, // save numpy matrix
        true  // save CSV
    )

    // Make file with smaple IDs and user-defined references or NA for each group
    samp_ref_pairs = sample_data
        .collectFile() { item -> [ "${item[4].id}.csv", "${item[0].id},${item[2].id ?: 'NA'}\n" ] }
        .map {[[id: it.getSimpleName()], it]}

    // For each group, assign references for variant calling if not user-defined
    ASSIGN_GROUP_REFERENCES (
        SOURMASH_COMPARE.out.csv.join(samp_ref_pairs).view()
    )

    emit:
    versions        = ch_versions                      // channel: [ versions.yml ]
}
