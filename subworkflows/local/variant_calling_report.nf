include { VCFTOTAB      } from '../../modules/local/vcftotab'

workflow VARIANT_CALLING_REPORT {

    take:
    ch_input // channel: [ val(ref_meta), vcf ]

    main:

    ch_versions = Channel.empty()

    VCFTOTAB ( ch_input )
    ch_versions = ch_versions.mix(VCFTOTAB.out.versions.first())
    VCFTOTAB.out.tab.view()

    emit:
    tab      = VCFTOTAB.out.tab           // channel: [ val(ref_meta), tab ]
    versions = ch_versions                     // channel: [ versions.yml ]
}

