include { PIRATE         } from '../../modules/nf-core/pirate/main'
include { SAMTOOLS_FAIDX } from '../../modules/nf-core/samtools/faidx/main'

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

    emit:
    pirate_aln      = PIRATE.out.aln        // channel: [ ref_meta, align_fasta ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

