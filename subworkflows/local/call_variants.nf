include { GRAPHTYPER_GENOTYPE       } from '../../modules/nf-core/graphtyper/genotype/main'
include { MAKEREGIONFILE            } from '../../modules/local/makeregionfile'
include { GRAPHTYPER_VCFCONCATENATE } from '../../modules/nf-core/graphtyper/vcfconcatenate/main'
include { TABIX_TABIX               } from '../../modules/nf-core/tabix/tabix/main'
include { GATK4_VARIANTFILTRATION   } from '../../modules/nf-core/gatk4/variantfiltration/main'
include { VCFLIB_VCFFILTER          } from '../../modules/nf-core/vcflib/vcffilter/main'

workflow CALL_VARIANTS {

    take:
    ch_input     // channel: [ val(meta), file(bam), file(bam_bai), val(ref_meta), file(reference), file(samtool_fai), file(picard_dict) ]

    main:

    ch_versions = Channel.empty()
    
    // group samples by reference genome
    //    ch_ref_grouped: [val(ref_meta), file(ref), file(samtools_fai), file(picard_dict), [val(meta)], [file(bam)],  [file(bam_bai)]]
    ch_ref_grouped = ch_input
        .groupTuple(by: 3)
        .map { [it[3], it[4][0], it[5][0], it[6][0], it[0], it[1], it[2]] } // remove redundant reference genome paths
    

    // make list of chromosome (fasta headers) names compatible with graphtyper
    ch_ref = ch_ref_grouped.map { it[0..1] }
    MAKEREGIONFILE ( ch_ref ) 

    // Run graphtyper on each group of samples for all chromosomes
    ch_ref_grouped = ch_ref_grouped.join(MAKEREGIONFILE.out.regions) // make inputs in same order
    GRAPHTYPER_GENOTYPE (
        ch_ref_grouped.map { [it[0], it[5], it[6]] },
        ch_ref_grouped.map { [it[0], it[1]] },
        ch_ref_grouped.map { [it[0], it[2]] },
        ch_ref_grouped.map { it[7] },
    )
    ch_versions = ch_versions.mix(GRAPHTYPER_GENOTYPE.out.versions.first())

    // Combine graphtyper VCFs for each group of samples
    GRAPHTYPER_VCFCONCATENATE ( GRAPHTYPER_GENOTYPE.out.vcf )
    ch_versions = ch_versions.mix(GRAPHTYPER_VCFCONCATENATE.out.versions.first())

    // Make tbi index for combined VCF
    TABIX_TABIX ( GRAPHTYPER_VCFCONCATENATE.out.vcf )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    // Filter heterozygous calls because bacteria are haploid, these are just errors
    vf_input = GRAPHTYPER_VCFCONCATENATE.out.vcf  // ensure inputs in same order
        .join(TABIX_TABIX.out.tbi)
        .join(ch_ref_grouped.map { it[0..3] })
    // vf_input: [val(ref_meta), file(vcf), file(tbi), file(ref), file(samtools_fai), file(picard_dict)]]
    GATK4_VARIANTFILTRATION (
        vf_input.map { it[0..2] },
        vf_input.map { it[3] },
        vf_input.map { it[4] },
        vf_input.map { it[5] }
    )
    ch_versions = ch_versions.mix(GATK4_VARIANTFILTRATION.out.versions.first())             

    // SelectVariants on the variant level, excluding non-variant sites:
    VCFLIB_VCFFILTER ( GATK4_VARIANTFILTRATION.out.vcf.join(GATK4_VARIANTFILTRATION.out.tbi) )
    ch_versions = ch_versions.mix(VCFLIB_VCFFILTER.out.versions.first()) 

    emit:
    vcf      = VCFLIB_VCFFILTER.out.vcf        // channel: [ ref_meta, vcf ]
    versions = ch_versions                     // channel: [ versions.yml ]
}

