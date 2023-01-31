include { TEST } from '../../modules/local/test'

workflow EUKARYOTEPIPELINE {

    take:
    input // channel: [ tuple val(meta), env(TAXON), path(hits), path(reads), path(reference) ]

    main:
    ch_taxon     = input.map { [it[0], it[1]] }
    ch_bbsketch  = input.map { [it[0], it[2]] }
    ch_reads     = input.map { [it[0], it[3]] }
    ch_reference = input.map { [it[0], it[4]] }
    ch_versions = Channel.empty()

    TEST( ch_taxon )

    emit:
    output = TEST.output

}
