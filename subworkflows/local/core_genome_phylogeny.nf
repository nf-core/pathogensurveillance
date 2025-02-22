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
include { BAKTA_BAKTA               } from '../../modules/nf-core/bakta/bakta/main'
include { BAKTA_BAKTADBDOWNLOAD     } from '../../modules/nf-core/bakta/baktadbdownload/main'

workflow CORE_GENOME_PHYLOGENY {

    take:
    sample_data
    ani_matrix // report_group_id, ani_matrix
    sample_assemblies // sample_id, assembly_path

    main:

    versions = Channel.empty()
    messages = Channel.empty()

    // Remove any samples that are not prokaryotes
    sample_data = sample_data
        .filter{it.kingdom == "Bacteria" || it.kingdom == "Archaea"}

    // Make file with sample IDs and user-defined references or NA for each group
    samp_ref_pairs = sample_data
        .map{ [it.sample_id, it.report_group_ids, it.ref_metas] }
        .transpose(by: 2)
        .map{ sample_id, report_group_id, ref_meta ->
            [sample_id, report_group_id, ref_meta.ref_id, ref_meta.ref_name, ref_meta.ref_description, ref_meta.ref_path, ref_meta.ref_contextual_usage]
        }
        .unique()
        .collectFile() { sample_id, report_group_id, ref_id, ref_name, ref_desc, ref_path, usage ->
            [ "${report_group_id}.tsv", "${sample_id}\t${ref_id}\t${ref_name}\t${ref_desc}\t${usage}\n" ]
        }
        .map {[[id: it.getSimpleName()], it]}

    // Assign referneces to groups for context in phylogenetic analyses
    ASSIGN_CONTEXT_REFERENCES (
        ani_matrix.combine(samp_ref_pairs, by: 0),
        params.n_ref_closest,
        params.n_ref_closest_named,
        params.n_ref_context
    )

    // Get relevant information from all references assigned to samples
    all_ref_data =  sample_data
        .map{ [[id: it.report_group_ids], it.ref_metas] }
        .transpose(by: 1)
        .map{ report_meta, ref_meta ->
            [[id: ref_meta.ref_id], report_meta, ref_meta.ref_path, ref_meta.gff]
        }
        .unique()

    // Get information for references selected for this analysis and check if they have an existing gff
    selected_ref_data = ASSIGN_CONTEXT_REFERENCES.out.references
        .splitText( elem: 1 )
        .map { [it[0], it[1].replace('\n', '')] } // remove newline that splitText adds
        .splitCsv( elem: 1, sep: '\t' )
        .map { report_meta, tsv_contents ->
            [[id: tsv_contents[0]], report_meta]
        }
        .join(all_ref_data, by: 0..1)
        .map { ref_meta, report_meta, ref_path, gff_path ->
            [ref_meta, report_meta, ref_path, gff_path]
        }
        .branch { ref_meta, report_meta, ref_path, gff_path ->
            has_gff: gff_path
            no_gff: ! gff_path
        }

    // Download the bakta database if needed
    //   Based on code from the bacass nf-core pipeline using the MIT license: https://github.com/nf-core/bacass
    if (params.bakta_db) {
        if (params.bakta_db.endsWith('.tar.gz')) {
            bakta_db_tar = Channel.fromPath(params.bakta_db).map{ [ [id: 'baktadb'], it] }
            UNTAR( bakta_db_tar )
            bakta_db = UNTAR.out.untar.map{ meta, db -> db }.first()
            versions = versions.mix(UNTAR.out.versions)
        } else {
            bakta_db = Channel.fromPath(params.bakta_db).first()
        }
    } else if (params.download_bakta_db) {
        BAKTA_BAKTADBDOWNLOAD()
        bakta_db = BAKTA_BAKTADBDOWNLOAD.out.db
        versions = versions.mix(BAKTA_BAKTADBDOWNLOAD.out.versions)
    }

    // Run bakta on samples and references without a gff already
    sample_assem_data = sample_data
        .map{ [[id: it.sample_id], [id: it.report_group_ids]] }
        .combine(sample_assemblies, by: 0)
    ref_assem_data = selected_ref_data.no_gff
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
        [] // prodigal_tf (optional)
    )
    versions = versions.mix(BAKTA_BAKTA.out.versions)

    // For references that have a gff already, combine the assembly with the gff the same way bakta outputs
    MAKE_GFF_WITH_FASTA (
        selected_ref_data.has_gff
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
    other_gffs = selected_ref_data.has_gff
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
    messages = messages.mix (
        SUBSET_CORE_GENES.out.message_data
            .splitCsv ( header:true, sep:'\t', quote:'"' )
            .map { file_name, message_data -> [
                message_data.sample_id ? [id: message_data.sample_id] : null,
                message_data.report_group_id ? [id: message_data.report_group_id] : null,
                message_data.reference_id ? [id: message_data.reference_id] : null,
                "CORE_GENOME_PHYLOGENY",
                message_data.message_type,
                message_data.description
            ] }
    )

    // Align each gene family with mafft
    core_genes = SUBSET_CORE_GENES.out.feat_seq // group_meta, [gene_dirs]
        .transpose() // group_meta, gene_dir
        .map { [[id: "${it[0].id}_${it[1].baseName}", group_id: it[0]], it[1]] } // subset_meta, gene_dir
    FILES_IN_DIR ( core_genes )
    MAFFT_SMALL ( FILES_IN_DIR.out.files.transpose(), [[], []], [[], []], [[], []], [[], []], [[], []], false )
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
    selected_refs = ASSIGN_CONTEXT_REFERENCES.out.references // group_meta, tsv
    versions      = versions             // versions.yml
    messages      = messages                // meta, group_meta, ref_meta, workflow, level, message

}
