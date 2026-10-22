include { SUBSET_SHARED_GENES } from '../../../modules/local/subset_shared_genes'
include { MAFFT_ALIGN         } from '../../../modules/nf-core/mafft/align'
include { IQTREE              } from '../../../modules/nf-core/iqtree'

workflow BUILD_GENE_TREES {

    take:
    gene_data_with_meta // group_meta, gene_seqs_dir, metadata

    main:

    messages = channel.empty()

    // Subset shared genes into clusters
    SUBSET_SHARED_GENES(gene_data_with_meta, params.phylo_min_genes, params.phylo_max_genes)

    // Construct warning messages from removed sample/reference IDs
    removed_with_meta = SUBSET_SHARED_GENES.out.removed_ids
        .join(gene_data_with_meta.map { group_meta, gene_seqs, metadata -> [group_meta, metadata] }, by: 0)

    messages = messages.mix(
        removed_with_meta
            .splitCsv(header: true, sep: '\t', quote: '"', elem: 1)
            .map { group_meta, row, metadata ->
                row.type == 'sample'
                    ? [[id: row.id], group_meta, null, metadata.workflow_name, 'WARNING', 'Sample removed from multigene phylogeny in order to find enough shared genes among other samples/references.']
                    : [null, group_meta, [id: row.id], metadata.workflow_name, 'NOTE', 'Reference removed from multigene phylogeny in order to find enough shared genes among other samples/references.']
            }
    )

    // Transpose to one tree per item and align
    tree_seqs = SUBSET_SHARED_GENES.out.feat_seqs
        .transpose()
        .map { group_meta, tree_seq_dir ->
            [
                [id: tree_seq_dir.baseName, group_id: group_meta],
                tree_seq_dir
            ]
        }
    MAFFT_ALIGN(tree_seqs, [[], []], [[], []], [[], []], [[], []], [[], []], false)

    // Infer phylogenetic tree from aligned shared genes
    iqtree_input = MAFFT_ALIGN.out.fas
        .map { meta, alignments ->
            [meta, alignments, []]
        }
    IQTREE(iqtree_input, [], [], [], [], [], [], [], [], [], [], [], [])

    // Group trees by report group
    trees = IQTREE.out.phylogeny
        .map { meta, tree -> [meta.group_id, tree] }
        .groupTuple(sort: 'hash')

    emit:
    phylogeny      = trees                                 // group_meta, [treefiles]
    gene_lists     = SUBSET_SHARED_GENES.out.gene_lists    // group_meta, [gene_list_tsvs]
    messages       = messages                              // meta, group_meta, ref_meta, workflow, level, message
    alignments     = MAFFT_ALIGN.out.fas                   // subset_meta, alignment

}
