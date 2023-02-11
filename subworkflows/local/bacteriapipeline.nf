
include { MAKE_REFERENCE_INDEX   } from './make_ref_index'  // this is being called from subworflows
include { ALIGN_READS_TO_REF     } from './align_reads_to_ref'
include { CALL_VARIANTS          } from './call_variants'
include { VARIANT_CALLING_REPORT } from './variant_calling_report'

workflow BACTERIAPIPELINE {

    take:
    input // channel: [ val(meta), val(taxon), path(hits), path(reads), val(ref_meta), path(reference) ]

    main:
    ch_taxon     = input.map { [it[0], it[1]] } 
    ch_bbsketch  = input.map { [it[0], it[2]] }
    ch_reads     = input.map { [it[0], it[3]] }
    ch_reference = input.map { [it[0], it[4], it[5]] }
    ch_versions = Channel.empty()

    MAKE_REFERENCE_INDEX ( ch_reference )
    
    ALIGN_READS_TO_REF (
        ch_reads
        .join(ch_reference)
        .join(MAKE_REFERENCE_INDEX.out.samtools_fai)
        .join(MAKE_REFERENCE_INDEX.out.bwa_index)
    )

    CALL_VARIANTS (
        ALIGN_READS_TO_REF.out.bam
        .join(ALIGN_READS_TO_REF.out.bai)
        .join(ch_reference)
        .join(MAKE_REFERENCE_INDEX.out.samtools_fai)
        .join(MAKE_REFERENCE_INDEX.out.picard_dict)
    )
    
    VARIANT_CALLING_REPORT (
        CALL_VARIANTS.out.vcf
    )

    emit:
    picard_dict  = MAKE_REFERENCE_INDEX.out.picard_dict
    samtools_fai = MAKE_REFERENCE_INDEX.out.samtools_fai
    samtools_gzi = MAKE_REFERENCE_INDEX.out.samtools_gzi
    bwa_index    = MAKE_REFERENCE_INDEX.out.bwa_index 
    versions = ch_versions                           // channel: [ versions.yml ]

}

