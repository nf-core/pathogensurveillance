include { GRAPHTYPER_GENOTYPE   } from '../../modules/nf-core/graphtyper/genotype/main'
include { MAKEREGIONFILE        } from '../../modules/local/makeregionfile'


workflow CALL_VARIANTS {

    take:
    ch_input     // channel: [ val(meta), file(bam), file(reference), val(ref_meta) ]

    main:

    ch_versions = Channel.empty()
    
    // group samples by reference genome
    //    ch_ref_grouped: [val(ref_meta), file(ref), [val(meta)], [file(bam)]]
    ch_ref_grouped = ch_input
        .groupTuple(by: 3)
        .map { [it[3], it[2][0], it[0], it[1]] } // remove redundant reference genome paths
    

    // make list of chromosome (fasta headers) names compatible with graphtyper
    ch_ref = ch_ref_grouped.map { it[0..1] }
    MAKEREGIONFILE ( ch_ref ) 
    MAKEREGIONFILE.out.regions.view()

    // Run graphtyper on each group of samples for all chromosomes

    // Combine graphtyper VCFs for each group of samples

    // Filter heterozygous calls because bacteria are haploid, these are just errors

    // SelectVariants on the variant level, excluding non-variant sites:

    emit:
    // vcf      = ?        // channel: [ reference, [ val(meta) ], vcf ]
    versions = ch_versions                     // channel: [ versions.yml ]
}

