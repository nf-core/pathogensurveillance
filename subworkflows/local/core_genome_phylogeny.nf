include { PIRATE                } from '../../modules/nf-core/pirate/main'
include { SAMTOOLS_FAIDX        } from '../../modules/nf-core/samtools/faidx/main'
include { MAFFT                 } from '../../modules/nf-core/mafft/main'
include { REFORMATPIRATERESULTS } from '../../modules/local/reformat_pirate_results'
include { ALIGNFEATURESEQUENCES } from '../../modules/local/align_feature_sequences'
include { SUBSETCOREGENES       } from '../../modules/local/subset_core_genes'

workflow CORE_GENOME_PHYLOGENY {

    take:
    ch_input // [ val(meta), file(gff), val(ref_meta), file(reference) ]

    main:

    ch_versions = Channel.empty()

    // group samples by reference genome                                        
    ch_gff_grouped = ch_input                                                   
        .groupTuple(by: 2)                                                      
        .map { [it[2], it[1]] } // remove redundant reference genome paths

    PIRATE ( ch_gff_grouped )
    ch_versions = ch_versions.mix(PIRATE.out.versions.first())

    REFORMATPIRATERESULTS ( PIRATE.out.results )                                                   
    ch_versions = ch_versions.mix(REFORMATPIRATERESULTS.out.versions.first())                  
    
    // Extract sequences of all genes (does not align, contrary to current name)
    ALIGNFEATURESEQUENCES ( PIRATE.out.results )                            
    ch_versions = ch_versions.mix(ALIGNFEATURESEQUENCES.out.versions.first())
    
    // Filter for core single copy genes and link their extracted sequences to a new folder
    SUBSETCOREGENES ( REFORMATPIRATERESULTS.out.gene_fam.join(ALIGNFEATURESEQUENCES.out.feat_seqs) )

    // Align each gene family with mafft
    MAFFT( SUBSETCOREGENES.out.feat_seq.transpose(), [] )                                         

    emit:
    pirate_aln      = PIRATE.out.aln        // channel: [ ref_meta, align_fasta ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

