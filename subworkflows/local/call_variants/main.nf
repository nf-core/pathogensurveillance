include { GRAPHTYPER_GENOTYPE       } from '../../../modules/local/graphtyper/genotype'
include { MAKE_REGION_FILE          } from '../../../modules/local/custom/make_region_file'
include { GRAPHTYPER_VCFCONCATENATE } from '../../../modules/local/graphtyper/vcfconcatenate'
include { TABIX_TABIX               } from '../../../modules/nf-core/tabix/tabix'
//include { TABIX_BGZIP               } from '../../../modules/local/bgzip/bgzip'
include { TABIX_BGZIP               } from '../../../modules/nf-core/tabix/bgzip'
include { GATK4_VARIANTFILTRATION   } from '../../../modules/nf-core/gatk4/variantfiltration'
include { VCFLIB_VCFFILTER          } from '../../../modules/nf-core/vcflib/vcffilter'

workflow CALL_VARIANTS {

    take:
    ch_input     // channel: [ val(meta), file(bam), file(bam_bai), val(ref_meta), file(reference), val(group_meta), file(samtool_fai), file(picard_dict) ]

    main:

    versions = Channel.empty()

    // group samples by reference genome and group
    //    ch_ref_grouped: [val(ref+group_meta), file(ref), file(samtools_fai), file(picard_dict), [val(meta)], [file(bam)],  [file(bam_bai)]]
    ch_ref_grouped = ch_input
        .map { [[id: "${it[5].id}_${it[3].id}", group: it[5], ref: it[3]], it[0], it[1], it[2], it[4], it[6], it[7]] }
        .groupTuple(by: 0, sort: 'hash')
        .map { [it[0], it[4].sort()[0], it[5].sort()[0], it[6].sort()[0], it[1], it[2], it[3]] } // remove redundant reference genome paths

    // make list of chromosome (fasta headers) names compatible with graphtyper
    ch_ref = ch_ref_grouped.map { it[0..1] }
    MAKE_REGION_FILE ( ch_ref )

    // Run graphtyper on each group of samples for all chromosomes
    ch_ref_grouped = ch_ref_grouped.join(MAKE_REGION_FILE.out.regions) // make inputs in same order
    GRAPHTYPER_GENOTYPE (
        ch_ref_grouped.map { [it[0], it[5], it[6]] },
        ch_ref_grouped.map { [it[0], it[1]] },
        ch_ref_grouped.map { [it[0], it[2]] },
        ch_ref_grouped.map { it[7] },
    )
    versions = versions.mix(GRAPHTYPER_GENOTYPE.out.versions)

    // Combine graphtyper VCFs for each group of samples
    GRAPHTYPER_VCFCONCATENATE ( GRAPHTYPER_GENOTYPE.out.vcf )
    versions = versions.mix(GRAPHTYPER_VCFCONCATENATE.out.versions)

    // Make tbi index for combined VCF
    TABIX_TABIX ( GRAPHTYPER_VCFCONCATENATE.out.vcf )
    versions = versions.mix(TABIX_TABIX.out.versions)

    // Make .gzi file from reference in case it is gzipped
    TABIX_BGZIP ( ch_ref )
    versions = versions.mix(TABIX_BGZIP.out.versions)

    // Filter heterozygous calls because bacteria are haploid, these are just errors
    vf_input = GRAPHTYPER_VCFCONCATENATE.out.vcf  //
        .join(TABIX_TABIX.out.tbi) // [val(ref+group_meta), file(vcf), file(tbi)]
        .join(ch_ref_grouped.map { it[0..3] })
        .join(TABIX_BGZIP.out.gzi) // [val(ref+group_meta), file(vcf), file(tbi), file(ref), file(samtools_fai), file(picard_dict), file(gzi)]]
    GATK4_VARIANTFILTRATION (
        vf_input.map { it[0..2] },
        vf_input.map { [it[0], it[3]] },
        vf_input.map { [it[0], it[4]] },
        vf_input.map { [it[0], it[5]] },
        vf_input.map { [it[0], it[6]] }
    )
    versions = versions.mix(GATK4_VARIANTFILTRATION.out.versions)

    // SelectVariants on the variant level, excluding non-variant sites:
    VCFLIB_VCFFILTER ( GATK4_VARIANTFILTRATION.out.vcf.join(GATK4_VARIANTFILTRATION.out.tbi) )
    versions = versions.mix(VCFLIB_VCFFILTER.out.versions)

    emit:
    vcf      = VCFLIB_VCFFILTER.out.vcf               // val(ref+group_meta), file(vcf)
    samples  = ch_ref_grouped.map { [it[0], it[4]] }  // val(ref+group_meta), [meta]
    versions = versions                               // versions.yml
}

