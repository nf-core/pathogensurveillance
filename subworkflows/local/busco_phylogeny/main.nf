include { BUSCO_BUSCO                                          } from '../../../modules/nf-core/busco/busco'
include { BUSCO_DOWNLOAD                                       } from '../../../modules/nf-core/busco/download'
include { ASSIGN_CONTEXT_REFERENCES as ASSIGN_BUSCO_REFERENCES } from '../../../modules/local/assign_context_references'
include { MAFFT_ALIGN as MAFFT_BUSCO                           } from '../../../modules/nf-core/mafft/align'
include { IQTREE as IQTREE_BUSCO                               } from '../../../modules/nf-core/iqtree'
include { SUBSET_BUSCO_GENES                                   } from '../../../modules/local/subset_busco_genes'

workflow BUSCO_PHYLOGENY {

    take:
    original_sample_data
    ani_matrix // report_group_id, ani_matrix
    sample_assemblies

    main:

    versions = Channel.empty()
    messages = Channel.empty()

    // Remove any samples that are not eukaryotes
    sample_data = original_sample_data
        .filter{it.domain == "Eukaryota"}

    // Make file with sample IDs and user-defined references or NA for each group
    samp_ref_pairs = sample_data
        .map{ [it.sample_id, it.report_group_ids, it.ref_metas] }
        .transpose(by: 2)
        .map{ sample_id, report_group_id, ref_meta ->
            [sample_id, report_group_id, ref_meta.ref_id, ref_meta.ref_name, ref_meta.ref_description, ref_meta.ref_path, ref_meta.ref_primary_usage]
        }
        .unique()
        .collectFile() { sample_id, report_group_id, ref_id, ref_name, ref_desc, ref_path, usage ->
            [ "${report_group_id}.tsv", "${sample_id}\t${ref_id}\t${ref_name}\t${ref_desc}\t${usage}\n" ]
        }
        .map {[[id: it.getSimpleName()], it]}

    // Assign referneces to groups for context in phylogenetic analyses
    ASSIGN_BUSCO_REFERENCES (
        ani_matrix.combine(samp_ref_pairs, by: 0),
        params.n_ref_closest,
        params.n_ref_closest_named,
        params.n_ref_context
    )
    versions = versions.mix(ASSIGN_BUSCO_REFERENCES.out.versions)

    // Create channel with required reference metadata and genomes from selected references
    references =  sample_data
        .map{ [[id: it.report_group_ids], it.ref_metas] }
        .transpose(by: 1)
        .map{ report_meta, ref_meta ->
            [report_meta, [id: ref_meta.ref_id], ref_meta.ref_path, ref_meta.ref_name]
        }
        .unique()

    selected_ref_data = ASSIGN_BUSCO_REFERENCES.out.references
        .splitText( elem: 1 )
        .map { [it[0], it[1].replace('\n', '')] } // remove newline that splitText adds
        .splitCsv( elem: 1, sep: '\t' )
        .map { report_meta, tsv_contents ->
            [report_meta, [id: tsv_contents[0]]]
        }
        .combine(references, by: 0..1)
        .map {report_meta, ref_meta, ref_path, ref_name ->
            [[id: ref_meta.id], report_meta, ref_path]
        }

    // Combine selected reference data with analogous sample metadata and assembled genomes
    busco_input = sample_data
        .map{ [[id: it.sample_id], [id: it.report_group_ids]] }
        .combine(sample_assemblies, by: 0)
        .mix(selected_ref_data)

    // Download BUSCO datasets
    BUSCO_DOWNLOAD ( Channel.from( "eukaryota_odb10" ) )
    versions = versions.mix(BUSCO_DOWNLOAD.out.versions)

    // Extract BUSCO genes for all unique reference genomes used in any sample/group
    BUSCO_BUSCO (
        busco_input
            .map{ meta, report_meta, path ->
                [meta, path]
            }
            .unique(),
        "genome",
        "eukaryota_odb10",
        BUSCO_DOWNLOAD.out.download_dir.first(), // .first() is needed to convert the queue channel to a value channel so it can be used multiple times.
        [],
        true
    )
    versions = versions.mix(BUSCO_BUSCO.out.versions)

    // Combine BUSCO output by gene for each report group
    sorted_busco_data = busco_input
        .combine(BUSCO_BUSCO.out.busco_dir, by: 0)
        .map{ meta, report_meta, path, busco_dir ->
            [report_meta, busco_dir]
        }
        .groupTuple(by: 0, sort: 'hash')
        .combine(samp_ref_pairs, by: 0)
    SUBSET_BUSCO_GENES (
        sorted_busco_data,
        params.phylo_min_genes,
        params.phylo_max_genes
    )
    versions = versions.mix(SUBSET_BUSCO_GENES.out.versions)
    messages = messages.mix (
        SUBSET_BUSCO_GENES.out.message_data
            .splitCsv ( header:true, sep:'\t', quote:'"' )
            .map { [
                it.sample_id == '' ? null : [id: it.sample_id],
                it.report_group_id == '' ? null : [id: it.report_group_id],
                it.reference_id == '' ? null : [id: it.reference_id],
                "BUSCO_PHYLOGENY",
                it.message_type,
                it.description
            ] }
    )

    // Align each gene family with mafft
    core_genes = SUBSET_BUSCO_GENES.out.feat_seqs
        .transpose()
        .map { report_meta, feat_seq_dir ->
            [
                [id: feat_seq_dir.baseName, group_id: report_meta],
                files(feat_seq_dir.resolve('*.*'), checkIfExists: true)
            ]
        }
        .transpose()
    MAFFT_BUSCO ( core_genes, [[], []], [[], []], [[], []], [[], []], [[], []], false )

    // Inferr phylogenetic tree from aligned core genes
    phylogeny_input = MAFFT_BUSCO.out.fas
        .groupTuple(sort: 'hash')
        .map { meta, alignments ->
            [meta, alignments, []]
        }
    IQTREE_BUSCO ( phylogeny_input, [], [], [], [], [], [], [], [], [], [], [], [] )
    versions = versions.mix(IQTREE_BUSCO.out.versions)
    trees = IQTREE_BUSCO.out.phylogeny // subset_meta, tree
        .map { [it[0].group_id, it[1]] } // group_meta, tree
        .groupTuple(sort: 'hash') // group_meta, [trees]

    emit:
    versions      = versions // versions.yml
    messages      = messages // meta, group_meta, ref_meta, workflow, level, message
    selected_refs = ASSIGN_BUSCO_REFERENCES.out.references
    tree          = trees

}
