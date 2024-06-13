include { KHMER_TRIMLOWABUND                          } from '../../modules/local/khmer_trimlowabund'
include { ASSIGN_MAPPING_REFERENCE                    } from '../../modules/local/assign_mapping_reference'
include { ASSIGN_CONTEXT_REFERENCES                   } from '../../modules/local/assign_context_references'
include { SOURMASH_SKETCH as SOURMASH_SKETCH_READS    } from '../../modules/nf-core/sourmash/sketch/main'
include { SOURMASH_SKETCH as SOURMASH_SKETCH_GENOME   } from '../../modules/nf-core/sourmash/sketch/main'
include { SOURMASH_COMPARE                            } from '../../modules/local/sourmash_compare'
include { SUBSET_READS                                } from '../../modules/local/subset_reads'

workflow SKETCH_COMPARISON {

    take:
    sample_data

    main:
    versions = Channel.empty()
    messages = Channel.empty()

    // Subset sample reads to increase speed of following steps
    SUBSET_READS (
        sample_data
            .map { [[id: it.sample_id[0]], it.paths[0], it.sendsketch_depth[0]] }
            .unique(),
        params.sketch_max_depth
    )
    versions = versions.mix(SUBSET_READS.out.versions.toSortedList().map{it[0]})

    // Trim rare k-mers from raw reads
    KHMER_TRIMLOWABUND (
        SUBSET_READS.out.reads
    )
    versions = versions.mix(KHMER_TRIMLOWABUND.out.versions)

    // Create signature for each sample
    SOURMASH_SKETCH_READS (
        KHMER_TRIMLOWABUND.out.sequence
    )
    versions = versions.mix(SOURMASH_SKETCH_READS.out.versions)

    // Create signature for each reference genome
    references = sample_data
        .map{ it.ref_metas }
        .transpose(by: 0)
        .map{ ref_meta -> [[id: ref_meta[0].ref_id], ref_meta[0].ref_path] }
        .unique()
    SOURMASH_SKETCH_GENOME (
        references
    )
    versions = versions.mix(SOURMASH_SKETCH_GENOME.out.versions)

    grouped_sample_sigs = sample_data
        .map { [[id: it.sample_id], [id: it.report_group_ids]] }
        .combine(SOURMASH_SKETCH_READS.out.signatures, by:0)
        .map { sample_id, report_group_id, signature -> [report_group_id, signature]}
        .unique()
    grouped_ref_sigs = sample_data
        .map{ [it.ref_metas, [id: it.report_group_ids]] }
        .transpose(by: 0)
        .map{ ref_meta, report_group_id -> [[id: ref_meta.ref_id], report_group_id] }
        .combine(SOURMASH_SKETCH_GENOME.out.signatures, by: 0)
        .map{ ref_id, report_group_id, signature -> [report_group_id, signature]}
        .unique()
    grouped_sigs = grouped_sample_sigs
        .mix(grouped_ref_sigs)
        .groupTuple()

    //// Make list of user-defined reference signatures for each group
    //user_sigs = sample_data
    //    .map { [it[2], it[4]] } // ref_meta, group_meta
    //    .combine(SOURMASH_SKETCH_GENOME.out.signatures, by: 0) // ref_meta, group_meta, sig
    //    .map { it[1..2] }  // [group_meta, sig], possibly duplicated
    //    .unique()
    //    .groupTuple() // group_meta, [sig]

    //// Make list of sample signatures for each group
    //sample_sigs = sample_data
    //    .combine(SOURMASH_SKETCH_READS.out.signatures, by: 0)
    //    .map { [it[4], it[5]] }
    //    .groupTuple() // group_meta, [sig]

    //// Make list of downloaded reference genome signatures for each group
    //assem_sigs = assem_samp_combos
    //    .combine(signatures, by: 0)
    //    .map { [it[1], it[0], it[2]] } // meta, assem, sig
    //    .combine(sample_data, by: 0) // meta, assem, sig, fastq, ref_meta, ref, group_meta
    //    .map { [it[6], it[2]] } // group_meta, sig
    //    .groupTuple() // group_meta, [sig]
    //    .map { [it[0], it[1].unique()] }

    //// Combine all signatures for each group
    //group_sigs = sample_sigs // group_meta, sample_sigs
    //    .join(assem_sigs) // group_meta, sample_sigs, assem_sigs
    //    .join(user_sigs, remainder: true) // group_meta, sample_sigs, assem_sigs, user_sigs
    //    .map { [it[0], it[1] + it[2] + (it[3] != null ? it[3] : [])] }

    // Compare all genomes/samples to eachother to create an ANI matrix
    SOURMASH_COMPARE (
        grouped_sigs,
        [], // file_list (optional)
        true, // save numpy matrix
        true  // save CSV
    )
    versions = versions.mix(SOURMASH_COMPARE.out.versions)

    //// Make file with smaple IDs and user-defined references or NA for each group
    //samp_ref_pairs = sample_data
    //    .collectFile() { item -> [ "${item[4].id}.csv", "${item[0].id},${item[2].id ?: 'NA'}\n" ] }
    //    .map {[[id: it.getSimpleName()], it]} // TODO this recreates the group_meta, but if other feilds besids "id" are added this will not preserve those

    //// For each group, assign references for variant calling if not user-defined
    //ASSIGN_MAPPING_REFERENCE (
    //    SOURMASH_COMPARE.out.csv.join(samp_ref_pairs),
    //    params.ref_min_ani
    //)

    //// Assign referneces to groups for context in phylogenetic analyses
    //ASSIGN_CONTEXT_REFERENCES (
    //    SOURMASH_COMPARE.out.csv.join(samp_ref_pairs),
    //    params.n_ref_closest,
    //    params.n_ref_context
    //)

    //// Convert CSV output back to nextflow channels
    //mapping_ref_ids = ASSIGN_MAPPING_REFERENCE.out.samp_ref_pairs
    //    .splitText( elem: 1 )
    //    .map { [it[0], it[1].replace('\n', '')] } // remove newline that splitText adds
    //    .splitCsv( elem: 1 )
    //    .map { [it[0].id] + it[1] } // [val(sample_id), val(group_id), val(reference_id)]
    //context_ref_ids = ASSIGN_CONTEXT_REFERENCES.out.samp_ref_pairs
    //    .splitText( elem: 1 )
    //    .map { [it[0], it[1].replace('\n', '')] } // remove newline that splitText adds
    //    .splitCsv( elem: 1 )
    //    .map { [it[0], it[1][0]] } // group_meta, ref_id

    //// Convert IDs back into full meta
    //id_meta_key = sample_data
    //    .map { [it[4].id, it[0].id, it[4], it[0]] }
    //mapping_refs = mapping_ref_ids
    //    .combine(id_meta_key, by: 0..1)
    //    .map { [it[4], it[3], it[2]] } // [val(meta), val(group_meta), val(ref_id)]
    //ref_id_key = assem_samp_combos // ref_meta, meta
    //    .map { ref_meta, meta -> [ref_meta.id, ref_meta] }

    //context_refs = context_ref_ids
    //    .map { group_meta, ref_id -> [ref_id, group_meta] }
    //    .combine(ref_id_key, by: 0)
    //    .map { ref_id, group_meta, ref_meta -> [group_meta, ref_meta] }
    //    .groupTuple()

    //// Add reference file based on ref_meta
    //user_refs = sample_data
    //    .filter { it[3] != null }
    //    .map { [it[2].id, it[2], it[3]] }
    //    .unique() // [val(ref_id), val(ref_meta), file(reference)]
    //null_refs = Channel.of ( ["__NULL__", [id: null], null])
    //all_refs = sequence
    //    .map { [it[0].id, it[0], it[1]] }
    //    .concat(user_refs) // [val(ref_id), val(ref_meta), file(reference)]
    //    .concat(null_refs)
    //mapping_refs_with_seq = mapping_refs
    //    .map { [it[2]] + it[0..1] } // [val(ref_id), val(meta), val(group_meta)]
    //    .combine(all_refs, by: 0)  // [val(ref_id), val(meta), val(group_meta), val(ref_meta), file(reference)]
    //    .map { it[1..4] }  // [val(meta), val(group_meta), val(ref_meta), file(reference)]

    //// Recreate sample data with new references picked
    //new_sample_data = sample_data  // meta, [reads], ref_meta, reference, group_meta
    //    .map { [it[0], it[4], it[1]] } // meta, group_meta, [reads]
    //    .combine ( assigned_refs_with_seq, by: 0..1 ) // meta, group_meta, [reads], ref_meta, reference
    //    .map { [it[0]] + it[2..4] + [it[1]] } // meta, [reads], ref_meta, reference, group_meta

    //// Report any samples that could not be assigned a reference
    //no_ref_warnings = new_sample_data // meta, [reads], ref_meta, reference, group_meta
    //    .filter { it[3] == null }
    //    .map { [it[0], it[4], null, "ASSIGN_REFERENCES", "WARNING", "Sample could not be assigned a reference, possibly because no similar orgnaism are present in NCBI RefSeq"] } // meta, group_meta, ref_meta, workflow, level, message
    //messages = messages.mix(no_ref_warnings)


    emit:
    ani_matrix    = SOURMASH_COMPARE.out.csv                   // group_meta, csv
    versions      = versions                                // versions
    messages      = messages                                   // meta, group_meta, ref_meta, workflow, level, message
}
