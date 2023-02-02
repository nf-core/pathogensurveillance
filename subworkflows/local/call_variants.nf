include { GRAPHTYPER_GENOTYPE                } from '../../modules/nf-core/graphtyper/genotype/main'


workflow CALL_VARIANTS {

    take:
    ch_input     // channel: [ val(meta), bam, reference ]

    main:

    ch_versions = Channel.empty()
    ch_input.view()
    
    // make list of chromosome (fasta headers) names compatible with graphtyper

    // group samples by reference genome

    // Run graphtyper on each group of samples for all chromosomes

    // Combine graphtyper VCFs for each group of samples

    // Filter heterozygous calls because bacteria are haploid, these are just errors

    // SelectVariants on the variant level, excluding non-variant sites:

    emit:
    // vcf      = ?        // channel: [ reference, [ val(meta) ], vcf ]
    versions = ch_versions                     // channel: [ versions.yml ]
}

