include { PIRATE                                              } from '../../../modules/nf-core/pirate'
include { MAFFT_ALIGN as MAFFT_CORE                           } from '../../../modules/nf-core/mafft/align'
include { IQTREE as IQTREE_CORE                               } from '../../../modules/nf-core/iqtree'
include { REFORMAT_PIRATE_RESULTS                             } from '../../../modules/local/reformat_pirate_results'
include { EXTRACT_FEATURE_SEQUENCES                           } from '../../../modules/local/extract_feature_sequences'
include { SUBSET_CORE_GENES                                   } from '../../../modules/local/subset_core_genes'
include { CALCULATE_POCP                                      } from '../../../modules/local/calculate_pocp'
include { ASSIGN_CONTEXT_REFERENCES as ASSIGN_CORE_REFERENCES } from '../../../modules/local/assign_context_references'
include { MAKE_GFF_WITH_FASTA                                 } from '../../../modules/local/make_gff_with_fasta'
include { BAKTA_BAKTA                                         } from '../../../modules/nf-core/bakta/bakta'
include { BAKTA_BAKTADBDOWNLOAD                               } from '../../../modules/nf-core/bakta/baktadbdownload'
include { UNTAR                                               } from '../../../modules/nf-core/untar'

workflow CORE_GENOME_PHYLOGENY {

    take:
    sample_data
    ani_matrix // report_group_id, ani_matrix
    sample_assemblies // sample_id, assembly_path

    main:

    messages = channel.empty()

    // Remove any samples that are not prokaryotes
    sample_data = sample_data
        .filter{ sample_meta -> sample_meta.domain == "Bacteria" || sample_meta.domain == "Archaea"}

    // Remove samples without a successful assembly
    sample_data = sample_data
        .map{ sample_meta -> [[id: sample_meta.sample_id], sample_meta] }
        .join(sample_assemblies, by: 0)
        .map{ sample_meta, sample_data_map, assembly_path -> sample_data_map }

    // Build stable raw TSV text of sample IDs and user-defined references for each group
    samp_ref_pairs = sample_data
        .map{ sample_meta -> [sample_meta.sample_id, sample_meta.report_group_ids, sample_meta.ref_metas] }
        .transpose(by: 2)
        .map{ sample_id, report_group_id, ref_meta ->
            [sample_id, report_group_id, ref_meta.ref_id, ref_meta.ref_name, ref_meta.ref_description, ref_meta.ref_path, ref_meta.ref_contextual_usage]
        }
        .unique()
        .map { sample_id, report_group_id, ref_id, ref_name, ref_desc, ref_path, usage ->
            [[id: report_group_id], "${sample_id}\t${ref_id}\t${ref_name}\t${ref_desc}\t${usage}"]
        }
        .groupTuple(by: 0, sort: 'hash')
        .map { group_meta, lines ->
            [group_meta, (lines.join('\n') + '\n').bytes.encodeBase64().toString()]
        }

    // Assign referneces to groups for context in phylogenetic analyses
    ASSIGN_CORE_REFERENCES (
        ani_matrix.combine(samp_ref_pairs, by: 0),
        params.n_ref_closest,
        params.n_ref_closest_named,
        params.n_ref_context
    )

    // Get relevant information from all references assigned to samples
    all_ref_data =  sample_data
        .map{ sample_meta -> [[id: sample_meta.report_group_ids], sample_meta.ref_metas] }
        .transpose(by: 1)
        .map{ report_meta, ref_meta ->
            [[id: ref_meta.ref_id], report_meta, ref_meta.ref_path, ref_meta.gff]
        }
        .unique()

    // Get information for references selected for this analysis and check if they have an existing gff
    selected_ref_data = ASSIGN_CORE_REFERENCES.out.references
        .splitText( elem: 1 )
        .map { row -> [row[0], row[1].replace('\n', '')] } // remove newline that splitText adds
        .splitCsv( elem: 1, sep: '\t' )
        .map { report_meta, tsv_contents ->
            [[id: tsv_contents[0].toString()], report_meta]
        }
        .join(all_ref_data, by: 0..1)
        .map { ref_meta, report_meta, ref_path, gff_path ->
            [ref_meta, report_meta, ref_path, gff_path]
        }
    selected_ref_data_has_gff = all_ref_data.filter { ref_meta, report_meta, ref_path, gff_path -> gff_path }
    selected_ref_data_no_gff = all_ref_data.filter { ref_meta, report_meta, ref_path, gff_path -> !gff_path }

    // Download the bakta database if needed
    //   Based on code from the bacass nf-core pipeline using the MIT license: https://github.com/nf-core/bacass
    if (params.bakta_db) {
        if (params.bakta_db.endsWith('.tar.gz')) {
            bakta_db_tar = channel.fromPath(params.bakta_db).map{ path -> [ [id: 'baktadb'], path] }
            UNTAR( bakta_db_tar )
            bakta_db = UNTAR.out.untar.map{ meta, db -> db }.first()
        } else {
            bakta_db = channel.fromPath(params.bakta_db).first()
        }
    } else if (params.download_bakta_db) {
        def bakta_cache = file("${params.data_dir}/bakta_db")
        if (bakta_cache.exists() && bakta_cache.listFiles()) {
            bakta_db = Channel.fromPath("${params.data_dir}/bakta_db/db-${params.bakta_db_type}").first()
        } else {
            BAKTA_BAKTADBDOWNLOAD()
            bakta_db = BAKTA_BAKTADBDOWNLOAD.out.db.first()
        }
    }

    // Run bakta on samples and references without a gff already
    sample_assem_data = sample_data
        .map{ sample_meta -> [[id: sample_meta.sample_id], [id: sample_meta.report_group_ids]] }
        .combine(sample_assemblies, by: 0)
    ref_assem_data = selected_ref_data_no_gff
        .map { ref_meta, report_meta, ref_path, gff_path ->
            [ref_meta, report_meta, ref_path]
        }
    all_assem_data = ref_assem_data
        .mix(sample_assem_data)
    BAKTA_BAKTA (
        all_assem_data
            .map { ref_meta, report_meta, assem_path ->
                [ref_meta, assem_path]
            }
            .unique(),
        bakta_db, // Bakta database
        [], // proteins (optional)
        [], // prodigal_tf (optional)
        [], // regions (optional)
        []  // hmms (optional)
    )

    // For references that have a gff already, combine the assembly with the gff the same way bakta outputs
    MAKE_GFF_WITH_FASTA (
        selected_ref_data_has_gff
            .map{ ref_meta, report_meta, ref_path, ref_gff ->
                [ref_meta, ref_path, ref_gff]
            }
    )

    // group samples by report group
    bakta_gffs = all_assem_data
        .combine(BAKTA_BAKTA.out.gff, by: 0) // sample_or_ref_meta, report_meta, assem_path, gff_path
        .map { sample_or_ref_meta, report_meta, assem_path, gff_path ->
            [sample_or_ref_meta, report_meta, gff_path]
        }
    other_gffs = selected_ref_data_has_gff
        .combine(MAKE_GFF_WITH_FASTA.out.gff, by: 0)
        .map{ ref_meta, report_meta, ref_path, ref_gff, ref_combined ->
            [ref_meta, report_meta, ref_combined]
        }
    all_gffs = bakta_gffs
        .mix(other_gffs)
        .map { sample_or_ref_meta, report_meta, gff ->
            [report_meta, gff]
        }
        .unique()
    PIRATE (
        all_gffs
            .groupTuple(by: 0, sort: 'hash')
            .map { meta, gffs -> // Attempt to fix intermittent "input file name collision" error
                def seen = [] as Set
                [meta, gffs.findAll { gff -> seen.add(gff.name) }]
            }
    )

    // Check that Pirate worked and report
    good_pirate_results = PIRATE.out.results
        .filter { result -> result[1].any{ path -> path.endsWith("PIRATE.gene_families.ordered.tsv") } }
    pirate_failed = PIRATE.out.results // val(group_meta), [result_files]
        .filter { result -> ! result[1].any{ path -> path.endsWith("PIRATE.gene_families.ordered.tsv") } }
        .map { result -> [null, result[0], null, "CORE_GENOME_PHYLOGENY", "WARNING", "Pirate failed to find a core genome, possibly becuase samples are very different or there are too few reads."] } // meta, group_meta, ref_meta, workflow, level, message
    messages = messages.mix(pirate_failed)

    REFORMAT_PIRATE_RESULTS ( good_pirate_results )


    // Calculate POCP from presence/absence matrix of genes
    CALCULATE_POCP (
        REFORMAT_PIRATE_RESULTS.out.gene_fam_pa
    )

    // Extract sequences of all genes
    EXTRACT_FEATURE_SEQUENCES ( good_pirate_results )

    // Filter for core single copy genes with no paralogs
    SUBSET_CORE_GENES (
        REFORMAT_PIRATE_RESULTS.out.gene_fam
            .join(EXTRACT_FEATURE_SEQUENCES.out.feat_seqs)
            .join(samp_ref_pairs),
        params.phylo_min_genes,
        params.phylo_max_genes
    )
    messages = messages.mix (
        SUBSET_CORE_GENES.out.message_data
            .splitCsv ( header:true, sep:'\t', quote:'"' )
            .map { file_name, message_data -> [
                message_data.sample_id ? [id: message_data.sample_id.toString()] : null,
                message_data.report_group_id ? [id: message_data.report_group_id.toString()] : null,
                message_data.reference_id ? [id: message_data.reference_id.toString()] : null,
                "CORE_GENOME_PHYLOGENY",
                message_data.message_type,
                message_data.description
            ] }
    )

    // Align each gene family with mafft
    core_genes = SUBSET_CORE_GENES.out.feat_seq // group_meta, [gene_dirs]
        .transpose() // group_meta, gene_dir
        .map { report_meta, feat_seq_dir ->
            [
                [id: feat_seq_dir.baseName, group_id: report_meta],
                feat_seq_dir
            ]
        }
    MAFFT_CORE ( core_genes, [[], []], [[], []], [[], []], [[], []], [[], []], false )

    // Inferr phylogenetic tree from aligned core genes
    phylogeny_input = MAFFT_CORE.out.fas
        .map { meta, alignments ->
            [meta, alignments, []]
        }
    IQTREE_CORE ( phylogeny_input, [], [], [], [], [], [], [], [], [], [], [], [] )
    trees = IQTREE_CORE.out.phylogeny // subset_meta, tree
        .map { meta -> [meta[0].group_id, meta[1]] } // group_meta, tree
        .groupTuple(sort: 'hash') // group_meta, [trees]

    // Mix in null placeholders for failed groups
    pirate_aln = pirate_failed // meta, group_meta, ref_meta, workflow, level, message
        .map { msg -> [msg[1], null] }
        .mix(PIRATE.out.aln) // group_meta, align_fasta
    phylogeny = pirate_failed // meta, group_meta, ref_meta, workflow, level, message
        .map { msg -> [msg[1], null] }
        .mix(trees) // group_meta, [trees]


    emit:
    pirate_aln    = pirate_aln              // group_meta, align_fasta
    phylogeny     = phylogeny               // group_meta, [trees]
    pocp          = CALCULATE_POCP.out.pocp // group_meta, pocp
    selected_refs = ASSIGN_CORE_REFERENCES.out.references // group_meta, tsv
    messages      = messages                // meta, group_meta, ref_meta, workflow, level, message

}
