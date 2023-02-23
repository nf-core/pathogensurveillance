include { VCFTOTAB      } from '../../modules/local/vcftotab'                   
include { VCFTOSNPALN   } from '../../modules/local/vcftosnpaln'                   
include { VARIANTREPORT } from '../../modules/local/variantreport'                

workflow VARIANT_CALLING_REPORT {

    take:
    ch_input // channel: [ val(ref_meta), vcf ]
    ch_samplesheet // channel: path

    main:

    ch_versions = Channel.empty()

    VCFTOTAB ( ch_input )
    ch_versions = ch_versions.mix(VCFTOTAB.out.versions.toSortedList().map{it[0]})

    VCFTOSNPALN ( VCFTOTAB.out.tab )
    ch_versions = ch_versions.mix(VCFTOSNPALN.out.versions.toSortedList().map{it[0]}) 

    VARIANTREPORT ( VCFTOSNPALN.out.fasta, ch_samplesheet )
    ch_versions = ch_versions.mix(VARIANTREPORT.out.versions.toSortedList().map{it[0]})             
               

    emit:
    tab      = VCFTOTAB.out.tab           // channel: [ val(ref_meta), tab ]
    versions = ch_versions                     // channel: [ versions.yml ]
}

