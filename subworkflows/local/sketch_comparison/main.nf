include { SOURMASH_SKETCH            } from '../../../modules/nf-core/sourmash/sketch/main'
include { SOURMASH_COMPARE           } from '../../../modules/nf-core/sourmash/compare'

workflow SKETCH_COMPARISON {

    take:
    sample_data
    assemblies

    main:
    versions = channel.empty()
    messages = channel.empty()

    // Create signature for each reference genome and successful assembly
    assemblies_to_sketch = sample_data
        .map{ [it.ref_metas] }
        .transpose(by: 0)
        .map{ ref_meta -> [[id: ref_meta[0].ref_id], ref_meta[0].ref_path] }
        .filter{ ref_id, ref_path -> ref_path }
        .unique()
        .mix(assemblies)
    SOURMASH_SKETCH (
        assemblies_to_sketch
    )
    versions = versions.mix(SOURMASH_SKETCH.out.versions)

    // Compare all genomes/samples to eachother to create an ANI matrix
    ref_sigs = sample_data
        .map{ [it.ref_metas, [id: it.report_group_ids]] }
        .transpose(by: 0)
        .map{ ref_meta, report_group_id -> [[id: ref_meta.ref_id], report_group_id] }
        .combine(SOURMASH_SKETCH.out.signatures, by: 0)
        .map{ ref_id, report_group_id, signature -> [report_group_id, signature]}
        .unique()
    assem_sigs = sample_data
        .map { [[id: it.sample_id], [id: it.report_group_ids]] }
        .combine(SOURMASH_SKETCH.out.signatures, by: 0)
        .map{ sample_id, report_group_id, signature -> [report_group_id, signature]}
        .unique()
    report_groups_with_assemblies = sample_data
        .map { [[id: it.sample_id], [id: it.report_group_ids]] }
        .join(assemblies.map { it[0] }, by: 0)
        .map { sample_id, report_group_id -> report_group_id }
        .unique()
    grouped_sigs = ref_sigs
        .mix(assem_sigs)
        .groupTuple(by: 0, sort: 'hash')
        .join(report_groups_with_assemblies, by: 0)
    SOURMASH_COMPARE (
        grouped_sigs,
        [], // file_list (optional)
        true, // save numpy matrix
        true  // save CSV
    )
    versions = versions.mix(SOURMASH_COMPARE.out.versions)

    emit:
    ani_matrix    = SOURMASH_COMPARE.out.csv                   // group_meta, csv
    versions      = versions                                // versions
    messages      = messages                                   // meta, group_meta, ref_meta, workflow, level, message
}
