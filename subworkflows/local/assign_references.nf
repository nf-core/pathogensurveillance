include { KHMER_TRIMLOWABUND                          } from '../../modules/local/khmer_trimlowabund'
include { ASSIGN_GROUP_REFERENCES                     } from '../../modules/local/assign_group_references'
include { SOURMASH_SKETCH as SOURMASH_SKETCH_READS    } from '../../modules/nf-core/sourmash/sketch/main'
include { SOURMASH_SKETCH as SOURMASH_SKETCH_GENOME   } from '../../modules/nf-core/sourmash/sketch/main'
include { SOURMASH_COMPARE                            } from '../../modules/local/sourmash_compare'
include { SUBSET_READS                                } from '../../modules/local/subset_reads'

workflow ASSIGN_REFERENCES {

    take:
    sample_data  // [val(meta), [file(fastq)], val(ref_meta), file(reference), val(group_meta)]
    assem_samp_combos // [val(assem_id), val(meta)] for each unique combination
    sequence // [val(assem_meta), file(fna)] for each assembly
    signatures // [val(assem_meta), file(sig)] for each assembly
    depth // [val(meta), val(depth)] for each sample

    main:
    ch_versions = Channel.empty()
    ch_messages = Channel.empty()

    // Subset sample reads to increase speed of following steps
    SUBSET_READS (
        sample_data                                                             
            .map { it[0..1] } // [val(meta), [file(fastq)]], possibly duplicated
            .combine(depth, by:0) // [val(meta), [file(fastq)], depth], possibly duplicated
            .unique(), // [val(meta), [file(fastq)], depth], one per sample                                                           
        params.sketch_max_depth                                                
    )                                                                           
    ch_versions = ch_versions.mix(SUBSET_READS.out.versions.toSortedList().map{it[0]})

    // Trim rare k-mers from raw reads
    KHMER_TRIMLOWABUND (
        SUBSET_READS.out.reads
    )
    ch_versions = ch_versions.mix(KHMER_TRIMLOWABUND.out.versions)

    // Create signature for each sample
    SOURMASH_SKETCH_READS (
        KHMER_TRIMLOWABUND.out.sequence
    )
    ch_versions = ch_versions.mix(SOURMASH_SKETCH_READS.out.versions)              
    
    // Create signature for each user-defined reference genome
    user_refs = sample_data
        .filter { it[3] != null }
        .map { it[2..3] } // val(ref_meta), file(ref)
        .unique()
    SOURMASH_SKETCH_GENOME (                                                           
        user_refs                                         
    )                                                                           
    ch_versions = ch_versions.mix(SOURMASH_SKETCH_GENOME.out.versions)           

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
        .join(user_sigs, remainder: true)
        .map { [it[0], it[1] + it[2] + (it[3] != null ? it[3] : [])] }

    // Compare all genomes/samples to eachother to create an ANI matrix
    SOURMASH_COMPARE (
        group_sigs,
        [], // file_list (optional)
        true, // save numpy matrix
        true  // save CSV
    )
    ch_versions = ch_versions.mix(SOURMASH_COMPARE.out.versions)          

    // Make file with smaple IDs and user-defined references or NA for each group
    samp_ref_pairs = sample_data
        .collectFile() { item -> [ "${item[4].id}.csv", "${item[0].id},${item[2].id ?: 'NA'}\n" ] }
        .map {[[id: it.getSimpleName()], it]} // TODO this recreates the group_meta, but if other feilds besids "id" are added this will not preserve those

    // For each group, assign references for variant calling if not user-defined
    ASSIGN_GROUP_REFERENCES (
        SOURMASH_COMPARE.out.csv.join(samp_ref_pairs)
    )

    // Convert CSV output back to nextflow channels
    assigned_refs_ids = ASSIGN_GROUP_REFERENCES.out.samp_ref_pairs
        .splitText( elem: 1 )
        .map { [it[0], it[1].replace('\n', '')] } // remove newline that splitText adds
        .splitCsv( elem: 1 )
        .map { [it[0].id] + it[1] } // [val(sample_id), val(group_id), val(reference_id)]
    

    // Convert IDs back into full meta
    id_meta_key = sample_data                                                                 
        .map { [it[4].id, it[0].id, it[4], it[0]] }
    assigned_refs = assigned_refs_ids
        .combine(id_meta_key, by: 0..1)
        .map { [it[4], it[3], it[2]] } // [val(meta), val(group_meta), val(ref_id)]

    // Add reference file based on ref_meta
    user_refs = sample_data
        .filter { it[3] != null }
        .map { [it[2].id, it[2], it[3]] } 
        .unique() // [val(ref_id), val(ref_meta), file(reference)]
    null_refs = Channel.of ( ["__NULL__", [id: null], null])
    all_refs = sequence
        .map { [it[0], [id: it[0]], it[1]] }
        .concat(user_refs) // [val(ref_id), val(ref_meta), file(reference)]
        .concat(null_refs)
    assigned_refs_with_seq = assigned_refs
        .map { [it[2]] + it[0..1] } // [val(ref_id), val(meta), val(group_meta)]
        .combine(all_refs, by: 0)  // [val(ref_id), val(meta), val(group_meta), val(ref_meta), file(reference)]
        .map { it[1..4] } // [val(meta), val(group_meta), val(ref_meta), file(reference)]

    // Recreate sample data with new references picked
    new_sample_data = sample_data
        .map { [it[0], it[4], it[1]] } // [val(meta), val(group_meta), [file(fastq)] ]
        .combine ( assigned_refs_with_seq, by: 0..1 ) // [val(meta), val(group_meta), [file(fastq)], val(ref_meta), file(reference)]
        .map { [it[0], it[2], it[3], it[4], it[1]] } // [val(meta), [file(fastq)], val(ref_meta), file(reference), val(group_meta)]


    emit:
    sample_data     = new_sample_data   // [val(meta), [file(fastq)], val(ref_meta), file(reference), val(group_meta)]
    ani_matrix      = SOURMASH_COMPARE.out.csv // [val(group_meta), val(csv)]
    versions        = ch_versions       // channel: [ versions.yml ]
}
