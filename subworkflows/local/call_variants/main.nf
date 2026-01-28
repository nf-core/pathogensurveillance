include { GRAPHTYPER_GENOTYPE       } from '../../../modules/nf-core/graphtyper/genotype'
include { GRAPHTYPER_VCFCONCATENATE } from '../../../modules/nf-core/graphtyper/vcfconcatenate'
include { TABIX_TABIX               } from '../../../modules/nf-core/tabix/tabix'
include { TABIX_BGZIP               } from '../../../modules/nf-core/tabix/bgzip'
include { GATK4_VARIANTFILTRATION   } from '../../../modules/nf-core/gatk4/variantfiltration'
include { VCFLIB_VCFFILTER          } from '../../../modules/nf-core/vcflib/vcffilter'

workflow CALL_VARIANTS {

    take:
    ch_input     // channel: [ val(meta), file(bam), file(bam_bai), val(ref_meta), file(reference), val(group_meta), file(samtool_fai), file(picard_dict) ]

    main:

    versions = Channel.empty()

    // group samples by reference genome and group
    ch_ref_grouped = ch_input
        .map { sample_meta, bam, bai, ref_meta, ref, report_meta, fai, dict ->
            [[id: "${report_meta.id}--${ref_meta.id}", group: report_meta, ref: ref_meta], sample_meta, bam, bai, ref, fai, dict]
        }
        .groupTuple(by: 0, sort: 'hash')
        .map { combined_meta, sample_meta, bam, bai, ref, fai, dict ->
            [combined_meta, sample_meta, bam, bai, ref.sort()[0], fai.sort()[0], dict.sort()[0]]
        }

    // make list of chromosome (fasta headers) names compatible with graphtyper
    combined_meta_key = ch_ref_grouped
        .map{ combined_meta, sample_meta, bam, bai, ref, fai, dict ->
            [combined_meta.id, combined_meta]
        }
        .unique()
    region_files = ch_ref_grouped
        .map { combined_meta, sample_meta, bam, bai, ref, fai, dict ->
            [combined_meta, ref]
        }
        .splitFasta( record: [id: true], elem: 1 )
        .collectFile(newLine: true) { combined_meta, header ->
            [ "${combined_meta.id}.txt", header.id ]
        }
        .map { [it.getSimpleName(), it] }
        .combine(combined_meta_key, by:0) // add back the full combined meta since collectFile only preserves the ID as a file name
        .map { combined_meta_id, region_path, combined_meta ->
            [combined_meta, region_path]
        }

    // Run graphtyper on each group of samples for all chromosomes
    ch_ref_grouped = ch_ref_grouped.join(region_files) // make inputs in same order
    GRAPHTYPER_GENOTYPE (
        ch_ref_grouped.map { combined_meta, sample_meta, bam, bai, ref, fai, dict, region_file ->
            [combined_meta, bam, bai]
        },
        ch_ref_grouped.map { combined_meta, sample_meta, bam, bai, ref, fai, dict, region_file ->
            [combined_meta, ref]
        },
        ch_ref_grouped.map { combined_meta, sample_meta, bam, bai, ref, fai, dict, region_file ->
            [combined_meta, fai]
        },
        ch_ref_grouped.map { combined_meta, sample_meta, bam, bai, ref, fai, dict, region_file ->
            [region_file]
        },
    )
    versions = versions.mix(GRAPHTYPER_GENOTYPE.out.versions)

    // Combine graphtyper VCFs for each group of samples
    GRAPHTYPER_VCFCONCATENATE ( GRAPHTYPER_GENOTYPE.out.vcf )
    versions = versions.mix(GRAPHTYPER_VCFCONCATENATE.out.versions)

    // Make tbi index for combined VCF
    TABIX_TABIX ( GRAPHTYPER_VCFCONCATENATE.out.vcf )

    // Make .gzi file from reference in case it is gzipped
    TABIX_BGZIP (
        ch_ref_grouped.map { combined_meta, sample_meta, bam, bai, ref, fai, dict, region_file ->
            [combined_meta, ref]
        }
    )
    versions = versions.mix(TABIX_BGZIP.out.versions)

    // Filter heterozygous calls because bacteria are haploid, these are just errors
    vf_input = GRAPHTYPER_VCFCONCATENATE.out.vcf  // make inputs in same order
        .join(TABIX_TABIX.out.index)
        .join(
            ch_ref_grouped.map { combined_meta, sample_meta, bam, bai, ref, fai, dict, region_file ->
                [combined_meta, ref, fai, dict]
            }
        )
        .join(TABIX_BGZIP.out.gzi) // [val(ref+group_meta), file(vcf), file(tbi), file(ref), file(fai), file(dict), file(gzi)]]
    GATK4_VARIANTFILTRATION (
        vf_input.map { combined_meta, vcf, tbi, ref, fai, dict, gzi ->
            [combined_meta, vcf, tbi]
        },
        vf_input.map { combined_meta, vcf, tbi, ref, fai, dict, gzi ->
            [combined_meta, ref]
        },
        vf_input.map { combined_meta, vcf, tbi, ref, fai, dict, gzi ->
            [combined_meta, fai]
        },
        vf_input.map { combined_meta, vcf, tbi, ref, fai, dict, gzi ->
            [combined_meta, dict]
        },
        vf_input.map { combined_meta, vcf, tbi, ref, fai, dict, gzi ->
            [combined_meta, gzi]
        }
    )
    versions = versions.mix(GATK4_VARIANTFILTRATION.out.versions)

    // SelectVariants on the variant level, excluding non-variant sites:
    VCFLIB_VCFFILTER ( GATK4_VARIANTFILTRATION.out.vcf.join(GATK4_VARIANTFILTRATION.out.tbi) )
    versions = versions.mix(VCFLIB_VCFFILTER.out.versions)

    emit:
    vcf      = VCFLIB_VCFFILTER.out.vcf
    samples  = ch_ref_grouped
        .map { combined_meta, sample_meta, bam, bai, ref, fai, dict, region_file ->
            [combined_meta, sample_meta]
        }
    versions = versions
}
