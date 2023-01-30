include { TEST } from '../../modules/local/test'

workflow EUKARYOTEPIPELINE {

    take:
    organism // channel: [ val(meta), [ bam ] ]

    main:
    ch_versions = Channel.empty()
    TEST ( organism )

    emit:
    output = TEST.output

}

