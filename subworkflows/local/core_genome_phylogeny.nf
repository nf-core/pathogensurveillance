include { PIRATE                    } from '../../modules/local/pirate.nf'
include { SAMTOOLS_FAIDX            } from '../../modules/nf-core/samtools/faidx/main'
include { MAFFT as MAFFT_SMALL      } from '../../modules/nf-core/mafft/main'
include { IQTREE2 as IQTREE2_CORE   } from '../../modules/local/iqtree2'
include { REFORMAT_PIRATE_RESULTS   } from '../../modules/local/reformat_pirate_results'
include { ALIGN_FEATURE_SEQUENCES   } from '../../modules/local/align_feature_sequences'
include { SUBSET_CORE_GENES         } from '../../modules/local/subset_core_genes'
include { RENAME_CORE_GENE_HEADERS  } from '../../modules/local/rename_core_gene_headers'
include { CALCULATE_POCP            } from '../../modules/local/calculate_pocp'

workflow CORE_GENOME_PHYLOGENY {

    take:
    sample_gff // [ val(meta), file(sample_gffs), val(group_meta), [file(ref_gffs)] ]
    ch_samplesheet // channel: path

    main:

    ch_versions = Channel.empty()
    messages = Channel.empty()

    // group samples by reference genome
    ch_gff_grouped = sample_gff
        .groupTuple(by: 2)
        .map { [it[2], it[1] + it[3].flatten().unique()] } // [ val(group_meta), [ gff ] ]

    PIRATE ( ch_gff_grouped )
    ch_versions = ch_versions.mix(PIRATE.out.versions.first())

    // Check that Pirate worked and report
    good_pirate_results = PIRATE.out.results
        .filter { it[1].any{ it.endsWith("PIRATE.gene_families.ordered.tsv") } }
    pirate_failed = PIRATE.out.results // val(group_meta), [result_files]
        .filter { ! it[1].any{ it.endsWith("PIRATE.gene_families.ordered.tsv") } }
        .map { [null, it[0], null, "CORE_GENOME_PHYLOGENY", "WARNING", "Pirate failed to find a core genome, possibly becuase samples are very different or there are too few reads."] } // meta, group_meta, ref_meta, workflow, level, message
    messages = messages.mix(pirate_failed)

    REFORMAT_PIRATE_RESULTS ( good_pirate_results )
    ch_versions = ch_versions.mix(REFORMAT_PIRATE_RESULTS.out.versions.first())


    // Calculate POCP from presence/absence matrix of genes
    CALCULATE_POCP (
       REFORMAT_PIRATE_RESULTS.out.gene_fam_pa
    )

    // Extract sequences of all genes (does not align, contrary to current name)
    ALIGN_FEATURE_SEQUENCES ( good_pirate_results )
    ch_versions = ch_versions.mix(ALIGN_FEATURE_SEQUENCES.out.versions.first())

    // Rename FASTA file headers to start with just sample ID for use with IQTREE
    RENAME_CORE_GENE_HEADERS ( ALIGN_FEATURE_SEQUENCES.out.feat_seqs )

    // Filter for core single copy genes with no paralogs
    SUBSET_CORE_GENES (
        REFORMAT_PIRATE_RESULTS.out.gene_fam.join(RENAME_CORE_GENE_HEADERS.out.feat_seqs),
        ch_samplesheet,
        params.min_core_genes,
        params.min_core_samps,
        params.min_core_refs,
        params.max_core_genes
    )

    // Align each gene family with mafft
    MAFFT_SMALL ( SUBSET_CORE_GENES.out.feat_seq.transpose(), [[], []], [[], []], [[], []], [[], []], [[], []] )
    ch_versions = ch_versions.mix(MAFFT_SMALL.out.versions.first())

    // Inferr phylogenetic tree from aligned core genes
    IQTREE2_CORE ( MAFFT_SMALL.out.fas.groupTuple(), [] )
    ch_versions = ch_versions.mix(IQTREE2_CORE.out.versions.first())

    // Mix in null placeholders for failed groups
    pirate_aln = pirate_failed // meta, group_meta, ref_meta, workflow, level, message
         .map { [it[1], null] }
         .mix(PIRATE.out.aln) // group_meta, align_fasta
    phylogeny = pirate_failed // meta, group_meta, ref_meta, workflow, level, message
         .map { [it[1], null] }
         .mix(IQTREE2_CORE.out.phylogeny) // group_meta, align_fasta


    emit:
    pirate_aln = pirate_aln              // group_meta, align_fasta
    phylogeny  = phylogeny               // group_meta, tree
    pocp       = CALCULATE_POCP.out.pocp // group_meta, pocp
    versions   = ch_versions             // versions.yml
    messages   = messages                // meta, group_meta, ref_meta, workflow, level, message

}
