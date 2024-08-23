include { BUSCO                     } from '../../modules/nf-core/busco/main'
include { BUSCO_DOWNLOAD            } from '../../modules/local/busco_download'
include { READ2TREE                 } from '../../modules/local/read2tree/main'
include { R2TF                      } from '../../modules/local/r2tf'
include { R2TDIR                    } from '../../modules/local/r2tdir'
include { R2TBIN                    } from '../../modules/local/r2tbin'
include { ASSIGN_CONTEXT_REFERENCES } from '../../modules/local/assign_context_references'
include { MAKE_READ2TREE_DB         } from '../../modules/local/make_read2tree_db'

workflow BUSCO_PHYLOGENY {

    take:
    original_sample_data
    ani_matrix // report_group_id, ani_matrix

    main:

    versions = Channel.empty()
    messages = Channel.empty()

    // Remove any samples that are not eukaryotes
    sample_data = original_sample_data
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
            [[id: ref_meta.id], report_meta, ref_path]
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
    grouped_busco_out = BUSCO.out.single_copy_fna
        .filter{ ref_meta, gene_paths ->
            gene_paths.size() > 0
        }
        .join(BUSCO.out.busco_dir)
        .combine(selected_ref_data, by: 0)
        .map { ref_meta, gene_paths, busco_dir, report_meta, ref_path ->
            [report_meta, busco_dir]
        }
        .unique()
        .groupTuple(by: 0)
    no_gene_warnings = BUSCO.out.single_copy_fna
        .filter{ ref_meta, gene_paths ->
            gene_paths.size() == 0
        }
        .combine(selected_ref_data, by: 0)
        .map { ref_meta, gene_paths, report_meta, ref_path ->
            [null, report_meta, ref_meta, "BUSCO_PHYLOGENY", "WARNING", "Reference excluded from BUSCO phylogeny because no single copy busco genes were found."]
        }
    messages = messages.mix(no_gene_warnings)
    MAKE_READ2TREE_DB ( grouped_busco_out, "eukaryota_odb10" )

    // group samples
    input_filtered = sample_data
        .map{ sample_meta ->
            [[id: sample_meta.sample_id], sample_meta.paths, [id: sample_meta.report_group_ids], sample_meta.sequence_type]
        }.unique()
    paired_end = input_filtered
        .filter { sample_meta, read_paths, report_meta, type ->
            (type == "illumina" || type == "bgiseq") && read_paths.size() == 2
        }
        .groupTuple(by: 2, sort: 'hash')
        .map { sample_metas, read_paths, report_meta, types ->
            [report_meta, sample_metas, read_paths.collect{it[0]}, read_paths.collect{it[1]}]
        }
    single_end = input_filtered
        .filter { sample_meta, read_paths, report_meta, type ->
            (type == "illumina" || type == "bgiseq") && read_paths.size() == 1
        }
        .groupTuple(by: 2, sort: 'hash')
        .map { sample_metas, read_paths, report_meta, types ->
            [report_meta, sample_metas, read_paths.collect{it[0]}]
        }
    longread = input_filtered
        .filter { sample_meta, read_paths, report_meta, type ->
            (type == "nanopore" || type == "pacbio") && read_paths.size() == 1
        }
        .groupTuple(by: 2, sort: 'hash')
        .map { sample_metas, read_paths, report_meta, types ->
            [report_meta, sample_metas, read_paths.collect{it[0]}]
        }

    // NOTE: This annoying way of adding parts of tuple channels one part at a time is becuase when a channel is empty,
    //   only one null is added, not the number of nulls equal to the the number of elements that should be in the tuple.
    r2t_input = MAKE_READ2TREE_DB.out.ref_aa
        .join(MAKE_READ2TREE_DB.out.ref_dna)
        .join(paired_end.map{report_meta, sample_metas, read_1, read_2 -> [report_meta, sample_metas]}, remainder:true)
        .join(paired_end.map{report_meta, sample_metas, read_1, read_2 -> [report_meta, read_1]}, remainder:true)
        .join(paired_end.map{report_meta, sample_metas, read_1, read_2 -> [report_meta, read_2]}, remainder:true)
        .join(single_end.map{report_meta, sample_metas, read_paths -> [report_meta, sample_metas]}, remainder:true)
        .join(single_end.map{report_meta, sample_metas, read_paths -> [report_meta, read_paths]}, remainder:true)
        .join(longread.map{report_meta, sample_metas, read_paths -> [report_meta, sample_metas]}, remainder:true)
        .join(longread.map{report_meta, sample_metas, read_paths -> [report_meta, read_paths]}, remainder:true)
        .map {report_meta, ref_aa, ref_dna, pair_meta, pair_1, pair_2, single_meta, single, long_meta, longreads ->
            [report_meta, pair_meta, pair_1, pair_2, single_meta, single, long_meta, longreads, ref_aa, ref_dna]
        }
        .map{ it.collect{ it ?: [] } } //replace nulls with empty lists

    // Run Read2tree
    READ2TREE ( r2t_input )

    emit:
    versions      = versions // versions.yml
    messages      = messages // meta, group_meta, ref_meta, workflow, level, message
    selected_refs = ASSIGN_CONTEXT_REFERENCES.out.references
    tree          = READ2TREE.out.tree // group_meta, tree
    r2t_ref_meta  = MAKE_READ2TREE_DB.out.ref_meta

}
