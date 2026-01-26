include { TRIM_AND_SKETCH            } from '../../../modules/local/trim_and_sketch'
include { SOURMASH_SKETCH            } from '../../../modules/nf-core/sourmash/sketch/main'
include { SOURMASH_COMPARE           } from '../../../modules/nf-core/sourmash/compare'

workflow SKETCH_COMPARISON {

    take:
    sample_data
    assemblies

    main:
    versions = Channel.empty()
    messages = Channel.empty()

    // Combine sample data with reads and assemblies to idenify which samples were not assembled
    sample_data
        .map { [[id: it.sample_id], it.paths] }
        .unique()
        .join(assemblies, remainder: true)
        .branch { sample_meta, read_paths, assembly ->
            reads: ! assembly
                return [sample_meta, read_paths]
            assemblies: assembly
                return [sample_meta, assembly]
        }
        .set { input }

    // Trim rare k-mers from raw reads and sketch
    TRIM_AND_SKETCH (
        input.reads
    )
    versions = versions.mix(TRIM_AND_SKETCH.out.versions)

    // Create signature for each reference genome
    assemblies_to_sketch = sample_data
        .map{ [it.ref_metas] }
        .transpose(by: 0)
        .map{ ref_meta -> [[id: ref_meta[0].ref_id], ref_meta[0].ref_path] }
        .unique()
        .mix(input.assemblies)
    SOURMASH_SKETCH (
        assemblies_to_sketch
    )
    versions = versions.mix(SOURMASH_SKETCH.out.versions)

    // Compare all genomes/samples to eachother to create an ANI matrix
    read_sigs = sample_data
        .map { [[id: it.sample_id], [id: it.report_group_ids]] }
        .combine(TRIM_AND_SKETCH.out.signatures, by:0)
        .map { sample_id, report_group_id, signature -> [report_group_id, signature]}
        .unique()
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
    grouped_sigs = read_sigs
        .mix(ref_sigs)
        .mix(assem_sigs)
        .groupTuple(by: 0, sort: 'hash')
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
