include { PIRATE                    } from '../../modules/local/pirate.nf'
include { SAMTOOLS_FAIDX            } from '../../modules/nf-core/samtools/faidx/main'
include { MAFFT as MAFFT_SMALL      } from '../../modules/nf-core/mafft/main'
include { IQTREE2 as IQTREE2_CORE   } from '../../modules/local/iqtree2'
include { REFORMAT_PIRATE_RESULTS   } from '../../modules/local/reformat_pirate_results'
include { ALIGN_FEATURE_SEQUENCES   } from '../../modules/local/align_feature_sequences'
include { SUBSET_CORE_GENES         } from '../../modules/local/subset_core_genes'
include { RENAME_CORE_GENE_HEADERS  } from '../../modules/local/rename_core_gene_headers'

workflow CORE_GENOME_PHYLOGENY {

    take:
    sample_gff // [ val(meta), file(sample_gffs), val(group_meta), [file(ref_gffs)] ]
    ch_samplesheet // channel: path

    main:

    ch_versions = Channel.empty()

    // group samples by reference genome                                        
    ch_gff_grouped = sample_gff                                                  
        .groupTuple(by: 2)                                                      
        .map { [it[2], it[1] + it[3].flatten().unique()] } // [ val(group_meta), [ gff ] ]

    PIRATE ( ch_gff_grouped )
    ch_versions = ch_versions.mix(PIRATE.out.versions.first())

    REFORMAT_PIRATE_RESULTS ( PIRATE.out.results )                                                   
    ch_versions = ch_versions.mix(REFORMAT_PIRATE_RESULTS.out.versions.first())                  
    
    // Extract sequences of all genes (does not align, contrary to current name)
    ALIGN_FEATURE_SEQUENCES ( PIRATE.out.results )                            
    ch_versions = ch_versions.mix(ALIGN_FEATURE_SEQUENCES.out.versions.first())

    // Rename FASTA file headers to start with just sample ID for use with IQTREE
    RENAME_CORE_GENE_HEADERS ( ALIGN_FEATURE_SEQUENCES.out.feat_seqs )
    
    // Filter for core single copy genes with no paralogs
    SUBSET_CORE_GENES ( REFORMAT_PIRATE_RESULTS.out.gene_fam.join(RENAME_CORE_GENE_HEADERS.out.feat_seqs) )

    // Align each gene family with mafft
    MAFFT_SMALL ( SUBSET_CORE_GENES.out.feat_seq.transpose(), [] )
    ch_versions = ch_versions.mix(MAFFT_SMALL.out.versions.first())

    // Inferr phylogenetic tree from aligned core genes
    IQTREE2_CORE ( MAFFT_SMALL.out.fas.groupTuple(), [] )
    ch_versions = ch_versions.mix(IQTREE2_CORE.out.versions.first())

    emit:
    pirate_aln      = PIRATE.out.aln             // channel: [ ref_meta, align_fasta ]
    phylogeny       = IQTREE2_CORE.out.phylogeny // channel: [ group_meta, tree ]
    versions        = ch_versions                // channel: [ versions.yml ]
}

