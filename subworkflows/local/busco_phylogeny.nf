include { BUSCO                     } from '../../modules/nf-core/busco/main'
include { BUSCO_DOWNLOAD            } from '../../modules/local/busco_download'
include { READ2TREE                 } from '../../modules/local/read2tree/main'
include { R2TF                      } from '../../modules/local/r2tf'
include { R2TDIR                    } from '../../modules/local/r2tdir'
include { R2TBIN                    } from '../../modules/local/r2tbin'
include { ASSIGN_CONTEXT_REFERENCES } from '../../modules/local/assign_context_references'

workflow BUSCO_PHYLOGENY {

    take:
    sample_data
    ani_matrix // report_group_id, ani_matrix
    // input // val(meta), [file(fastq)], val(group_meta), [val(ref_meta)], val(kingdom), val(depth)
    //ref_seqs // val(ref_meta), file(fna)

    main:

    versions = Channel.empty()
    messages = Channel.empty()

    // Remove any samples that are not eukaryotes
    sample_data = sample_data
        .filter{it.kingdom == "Eukaryota"}

    // Make file with sample IDs and user-defined references or NA for each group
    samp_ref_pairs = sample_data
        .map{ [[id: it.sample_id], [id: it.report_group_ids], it.ref_metas] }
        .transpose(by: 2)
        .map{ sample_meta, report_meta, ref_meta ->
            [sample_meta, report_meta, [id: ref_meta.ref_id], ref_meta.ref_path, ref_meta.ref_primary_usage]
        }
        .collectFile() { sample_meta, report_meta, ref_id, ref_path, usage ->
            [ "${report_meta.id}.csv", "${sample_meta.id},${ref_id.id},${usage}\n" ]
        }
        .map {[[id: it.getSimpleName()], it]}

    // Assign referneces to groups for context in phylogenetic analyses
    ASSIGN_CONTEXT_REFERENCES (
        ani_matrix.combine(samp_ref_pairs, by: 0),
        params.n_ref_closest,
        params.n_ref_context
    )
    references =  sample_data
        .map{ [[id: it.report_group_ids], it.ref_metas] }
        .transpose(by: 1)
        .map{ report_meta, ref_meta ->
            [report_meta, [id: ref_meta.ref_id], ref_meta.ref_path, ref_meta.ref_name]
        }
    selected_ref_data = ASSIGN_CONTEXT_REFERENCES.out.references
        .splitText( elem: 1 )
        .map { [it[0], it[1].replace('\n', '')] } // remove newline that splitText adds
        .splitCsv( elem: 1 )
        .map { report_meta, csv_contents ->
            [report_meta, [id: csv_contents[0]]]
        }
        .combine(references, by: 0..1)
        .map {report_meta, ref_meta, ref_path, ref_name ->
            [[id: ref_meta.id, name: ref_name], report_meta, ref_path]
        }

    // Only use reference with unique species names TODO: this relies on users to name reference correctly. A better solution should be found
    selected_ref_data = selected_ref_data
        .unique { ref_meta, report_meta, ref_path ->
            ["${ref_meta.name.substring(0,3)}${ref_meta.name.split(' ')[1].substring(0, 2)}".toUpperCase(), report_meta]
        }

    // Download BUSCO datasets
    BUSCO_DOWNLOAD ( Channel.from( "eukaryota_odb10" ) )

    // Extract BUSCO genes for all unique reference genomes used in any sample/group
    BUSCO (
        selected_ref_data
            .map{ ref_meta, report_meta, ref_path ->
                [ref_meta, ref_path]
            }
            .unique(),
        "genome",
        "eukaryota_odb10",
        BUSCO_DOWNLOAD.out.download_dir.first(), // .first() is needed to convert the queue channel to a value channel so it can be used multiple times.
        []
    )

    // Create Read2tree database
    R2TF (
        BUSCO.out.busco_dir
	)

    // Create directories for Read2Tree
    R2TDIR (
        R2TF.out.output
            .combine(selected_ref_data, by: 0)
            .map { ref_meta, rt2f_out, report_meta, ref_path ->
                [report_meta, rt2f_out]
            }
            .groupTuple(by: 0, sort: 'hash')
    )

    // Bin busco files by genes
    R2TBIN (
        R2TDIR.out.markers
    )

    // group samples
    input_filtered = sample_data
        .map{ sample_meta ->
            [[id: sample_meta.sample_id], sample_meta.paths, [id: sample_meta.report_group_ids], sample_meta.sequence_type]
        }
        .unique()
    paired_end = input_filtered
        .filter { sample_meta, read_paths, report_meta, type ->
            type == "illumina" && read_paths.size() == 2
        }
        .groupTuple(by: 2)
        .map { sample_metas, read_paths, report_meta, types ->
            [report_meta, read_paths.collect{it[0]}, read_paths.collect{it[1]}]
        }
    single_end = input_filtered
        .filter { sample_meta, read_paths, report_meta, type ->
            type == "illumina" && read_paths.size() == 1
        }
        .groupTuple(by: 2)
        .map { sample_metas, read_paths, report_meta, types ->
            [report_meta, reads_paths]
        }
    longread = input_filtered
        .filter { sample_meta, read_paths, report_meta, type ->
            type == "nanopore" && read_paths.size() == 1
        }
        .groupTuple(by: 2)
        .map { sample_metas, read_paths, report_meta, types ->
            [report_meta, reads_paths]
        }
    rt2_input = paired_end
        .join(single_end, remainder:true)// group_meta, pair_1, pair_2, single_end
        .join(longread, remainder:true) // group_meta, pair_1, pair_2, single_end, longread
        .map { [it[0], it[1] ?: [], it[2] ?: [], it[3] ?: [], it[4] ?: []] } // replace nulls with []
        .join( R2TBIN.out.busco_markers )
        .join( R2TDIR.out.dna_ref )

    // Run Read2tree
    READ2TREE ( rt2_input )

    emit:
    versions      = versions // versions.yml
    messages      = messages // meta, group_meta, ref_meta, workflow, level, message
    selected_refs = ASSIGN_CONTEXT_REFERENCES.out.references
    tree          = READ2TREE.out.tree // group_meta, tree

}
