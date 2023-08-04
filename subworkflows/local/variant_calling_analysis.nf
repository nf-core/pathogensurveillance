include { MAKE_REFERENCE_INDEX   } from './make_ref_index'  // this is being called from subworkflows
include { ALIGN_READS_TO_REF     } from './align_reads_to_ref'
include { CALL_VARIANTS          } from './call_variants'
include { VARIANT_CALLING_REPORT } from './variant_calling_report'

workflow VARIANT_CALLING_ANALYSIS {

    take:
    input // [val(meta), [file(fastq)], val(ref_meta), file(reference), val(group_meta)]
    ch_samplesheet // channel: path

    main:
    ch_reads     = input.map { [it[0], it[1]] }
    ch_reference = input.map { [it[0], it[2], it[3]] }
    ch_versions = Channel.empty()

    MAKE_REFERENCE_INDEX ( ch_reference.unique() )
    ch_versions = ch_versions.mix(MAKE_REFERENCE_INDEX.out.versions)      
    
    ALIGN_READS_TO_REF (
        ch_reads
        .join(ch_reference)
        .join(MAKE_REFERENCE_INDEX.out.samtools_fai)
        .join(MAKE_REFERENCE_INDEX.out.bwa_index)
        .unique()
    )
    ch_versions = ch_versions.mix(ALIGN_READS_TO_REF.out.versions)       

    CALL_VARIANTS (
        ALIGN_READS_TO_REF.out.bam
        .join(ALIGN_READS_TO_REF.out.bai, by: 0..1) // [meta, ref_meta, bam, bai]
        .combine(input.map { [it[0], it[2], it[3], it[4]] }, by: 0..1)  // [meta, ref_meta, bam, bai, ref, group_meta]
        .combine(MAKE_REFERENCE_INDEX.out.samtools_fai, by: 1)
        .combine(MAKE_REFERENCE_INDEX.out.picard_dict, by: 1)
    )
    ch_versions = ch_versions.mix(CALL_VARIANTS.out.versions)       
    
    VARIANT_CALLING_REPORT (
        CALL_VARIANTS.out.vcf,
        ch_samplesheet
    )
    ch_versions = ch_versions.mix(VARIANT_CALLING_REPORT.out.versions)       

    emit:
    picard_dict  = MAKE_REFERENCE_INDEX.out.picard_dict
    samtools_fai = MAKE_REFERENCE_INDEX.out.samtools_fai
    samtools_gzi = MAKE_REFERENCE_INDEX.out.samtools_gzi
    bwa_index    = MAKE_REFERENCE_INDEX.out.bwa_index 
    versions = ch_versions                           // channel: [ versions.yml ]

}

