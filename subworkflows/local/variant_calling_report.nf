include { VCFTOTAB      } from '../../modules/local/vcftotab'                   
include { VCFTOSNPALN   } from '../../modules/local/vcftosnpaln'                   
include { VARIANTREPORT } from '../../modules/local/variantreport'                

workflow VARIANT_CALLING_REPORT {

    take:
    ch_input // channel: [ val(ref_meta), vcf ]

    main:

    ch_versions = Channel.empty()

    VCFTOTAB ( ch_input )
    ch_versions = ch_versions.mix(VCFTOTAB.out.versions.first())

    VCFTOSNPALN ( VCFTOTAB.out.tab )
    ch_versions = ch_versions.mix(VCFTOSNPALN.out.versions.first()) 

    VARIANTREPORT ( VCFTOSNPALN.out.fasta )
    ch_versions = ch_versions.mix(VARIANTREPORT.out.versions.first())             
               

    emit:
    tab      = VCFTOTAB.out.tab           // channel: [ val(ref_meta), tab ]
    versions = ch_versions                     // channel: [ versions.yml ]
}

