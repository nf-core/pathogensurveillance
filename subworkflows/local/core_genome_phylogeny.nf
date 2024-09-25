include { PIRATE                    } from '../../modules/local/pirate.nf'
include { SAMTOOLS_FAIDX            } from '../../modules/nf-core/samtools/faidx/main'
include { MAFFT as MAFFT_SMALL      } from '../../modules/nf-core/mafft/main'
include { IQTREE2 as IQTREE2_CORE   } from '../../modules/local/iqtree2'
include { REFORMAT_PIRATE_RESULTS   } from '../../modules/local/reformat_pirate_results'
include { ALIGN_FEATURE_SEQUENCES   } from '../../modules/local/align_feature_sequences'
include { SUBSET_CORE_GENES         } from '../../modules/local/subset_core_genes'
include { RENAME_CORE_GENE_HEADERS  } from '../../modules/local/rename_core_gene_headers'
include { CALCULATE_POCP            } from '../../modules/local/calculate_pocp'
include { FILES_IN_DIR              } from '../../modules/local/files_in_dir.nf'
include { ASSIGN_CONTEXT_REFERENCES } from '../../modules/local/assign_context_references'
include { MAKE_GFF_WITH_FASTA       } from '../../modules/local/make_gff_with_fasta.nf'

workflow CORE_GENOME_PHYLOGENY {

    take:
    sample_data
    ani_matrix // report_group_id, ani_matrix
    sample_and_ref_gff // sample_id, gff

    main:

    versions = Channel.empty()
    messages = Channel.empty()

    // Remove any samples that are not prokaryotes
    sample_data = sample_data
        .filter{it.kingdom == "Bacteria" || it.kingdom == "Archaea"}

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
            [[id: ref_meta.ref_id], report_meta, ref_meta.ref_path, ref_meta.gff]
        }
        .combine(sample_and_ref_gff, by: 0)
        .map { ref_id, report_id, ref_path, downloaded_gff, new_gff ->
            [report_id, ref_id, ref_path, downloaded_gff ?: new_gff]
        }

    // Combine refernece sequence with reference gffs for use with pirate
    selected_ref_data = ASSIGN_CONTEXT_REFERENCES.out.references
        .splitText( elem: 1 )
        .map { [it[0], it[1].replace('\n', '')] } // remove newline that splitText adds
        .splitCsv( elem: 1 )
        .map { report_meta, csv_contents ->
            [report_meta, [id: csv_contents[0]]]
        }
        .join(references, by: 0..1)
        .map {report_meta, ref_meta, ref_path, ref_gff ->
            [ref_meta, report_meta, ref_path, ref_gff]
        }
    MAKE_GFF_WITH_FASTA (
        selected_ref_data
            .map{ ref_meta, report_meta, ref_path, ref_gff ->
                [ref_meta, ref_path, ref_gff]
            }
            .unique()
    )

    // group samples by report group
    ref_gff_data = selected_ref_data
        .combine(MAKE_GFF_WITH_FASTA.out.gff, by: 0)
        .map{ ref_meta, report_meta, ref_path, ref_gff, ref_combined ->
            [ref_meta, report_meta, ref_combined]
        }
    sample_and_ref_gff_data = sample_data
        .map{ [[id: it.sample_id], [id: it.report_group_ids]] }
        .combine(sample_and_ref_gff, by: 0) // sample_meta, report_meta, gff
    gff_data = sample_and_ref_gff_data
        .mix(ref_gff_data)
        .map { sample_meta, report_meta, gff ->
            [report_meta, gff]
        }
        .unique()
    PIRATE (
        gff_data
            .groupTuple(by: 0, sort: 'hash')
    )
    versions = versions.mix(PIRATE.out.versions)

    // Check that Pirate worked and report
    good_pirate_results = PIRATE.out.results
        .filter { it[1].any{ it.endsWith("PIRATE.gene_families.ordered.tsv") } }
    pirate_failed = PIRATE.out.results // val(group_meta), [result_files]
        .filter { ! it[1].any{ it.endsWith("PIRATE.gene_families.ordered.tsv") } }
        .map { [null, it[0], null, "CORE_GENOME_PHYLOGENY", "WARNING", "Pirate failed to find a core genome, possibly becuase samples are very different or there are too few reads."] } // meta, group_meta, ref_meta, workflow, level, message
    messages = messages.mix(pirate_failed)

    REFORMAT_PIRATE_RESULTS ( good_pirate_results )
    versions = versions.mix(REFORMAT_PIRATE_RESULTS.out.versions)


    // Calculate POCP from presence/absence matrix of genes
    CALCULATE_POCP (
       REFORMAT_PIRATE_RESULTS.out.gene_fam_pa
    )

    // Extract sequences of all genes (does not align, contrary to current name)
    ALIGN_FEATURE_SEQUENCES ( good_pirate_results )
    versions = versions.mix(ALIGN_FEATURE_SEQUENCES.out.versions)

    // Rename FASTA file headers to start with just sample ID for use with IQTREE
    RENAME_CORE_GENE_HEADERS ( ALIGN_FEATURE_SEQUENCES.out.feat_seqs )

    // Filter for core single copy genes with no paralogs
    SUBSET_CORE_GENES (
        REFORMAT_PIRATE_RESULTS.out.gene_fam
            .join(RENAME_CORE_GENE_HEADERS.out.feat_seqs)
            .join(samp_ref_pairs),
        params.phylo_min_genes,
        params.phylo_max_genes
    )

    // Report any sample or references that have been removed from the analysis
    removed_refs = SUBSET_CORE_GENES.out.removed_ref_ids
        .splitText()
        .map { [null, [id: it[1].replace('\n', '')], it[0], "CORE_GENOME_PHYLOGENY", "WARNING", "Reference removed from core gene phylogeny in order to find enough core genes."] } // meta, group_meta, ref_meta, workflow, level, message
    removed_samps = SUBSET_CORE_GENES.out.removed_sample_ids
        .splitText()
        .map { [[id: it[1].replace('\n', '')], null, it[0], "CORE_GENOME_PHYLOGENY", "WARNING", "Sample removed from core gene phylogeny in order to find enough core genes."] } // meta, group_meta, ref_meta, workflow, level, message
    messages = messages.mix(removed_refs)
    messages = messages.mix(removed_samps)

    // Align each gene family with mafft
    core_genes = SUBSET_CORE_GENES.out.feat_seq // group_meta, [gene_dirs]
        .transpose() // group_meta, gene_dir
        .map { [[id: "${it[0].id}_${it[1].baseName}", group_id: it[0]], it[1]] } // subset_meta, gene_dir
    FILES_IN_DIR ( core_genes )
    MAFFT_SMALL ( FILES_IN_DIR.out.files.transpose(), [[], []], [[], []], [[], []], [[], []], [[], []] )
    versions = versions.mix(MAFFT_SMALL.out.versions)

    // Inferr phylogenetic tree from aligned core genes
    IQTREE2_CORE ( MAFFT_SMALL.out.fas.groupTuple(sort: 'hash'), [] )
    versions = versions.mix(IQTREE2_CORE.out.versions)
    trees = IQTREE2_CORE.out.phylogeny // subset_meta, tree
        .map { [it[0].group_id, it[1]] } // group_meta, tree
        .groupTuple(sort: 'hash') // group_meta, [trees]

    // Mix in null placeholders for failed groups
    pirate_aln = pirate_failed // meta, group_meta, ref_meta, workflow, level, message
         .map { [it[1], null] }
         .mix(PIRATE.out.aln) // group_meta, align_fasta
    phylogeny = pirate_failed // meta, group_meta, ref_meta, workflow, level, message
         .map { [it[1], null] }
         .mix(trees) // group_meta, [trees]


    emit:
    pirate_aln    = pirate_aln              // group_meta, align_fasta
    phylogeny     = phylogeny               // group_meta, [trees]
    pocp          = CALCULATE_POCP.out.pocp // group_meta, pocp
    selected_refs = ASSIGN_CONTEXT_REFERENCES.out.references // group_meta, csv
    versions      = versions             // versions.yml
    messages      = messages                // meta, group_meta, ref_meta, workflow, level, message

}
