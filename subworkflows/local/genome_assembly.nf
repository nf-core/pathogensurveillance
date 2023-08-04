include { FASTP           } from '../../modules/nf-core/fastp/main'
include { SPADES          } from '../../modules/nf-core/spades/main'
include { FILTER_ASSEMBLY } from '../../modules/local/filter_assembly'
include { QUAST           } from '../../modules/nf-core/quast/main'
include { BAKTA_BAKTA     } from '../../modules/nf-core/bakta/bakta/main'

workflow GENOME_ASSEMBLY {

    take:
    ch_input // channel: [ val(meta), [fastq_1, fastq_2], val(ref_meta), file(reference), val(group_meta), val(kingdom) ]

    main:

    ch_versions = Channel.empty()
    ch_input_filtered = ch_input
        .filter { it[5] == "Bacteria" }
    ch_reads = ch_input_filtered
        .map { it[0..1] }
        .unique()
    
    FASTP ( ch_reads, [], false, false )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    SPADES (
        FASTP.out.reads.map { [it[0], it[1], [], []] }, // tuple val(meta), path(illumina), path(pacbio), path(nanopore)
        [], // val yml
        []  // val hmm
    )
    ch_versions = ch_versions.mix(SPADES.out.versions.first())

    FILTER_ASSEMBLY (
        SPADES.out.scaffolds
    )
    ch_versions = ch_versions.mix(FILTER_ASSEMBLY.out.versions.first())

    ch_ref_grouped = ch_input_filtered
        .combine(FILTER_ASSEMBLY.out.filtered)
        .groupTuple(by: 2)
        .map { [it[2], it[3].sort()[0], it[1]] }
    QUAST (
        ch_ref_grouped.map { it[2] }, // consensus (one or more assemblies)
        ch_ref_grouped.map { it[1] }, // fasta (reference, optional)
        [], // gff (optional)
        true, // use_fasta
        false // use_gff
    )

    ch_bakta_db = Channel.value("$projectDir/assets/bakta_db/db-light")
    BAKTA_BAKTA (
        FILTER_ASSEMBLY.out.filtered, // Genome assembly
        ch_bakta_db, // Bakta database
        [], // proteins (optional)
        [] // prodigal_tf (optional)
    )

    emit:
    reads     = FASTP.out.reads           // channel: [ val(meta), [ fastq_1, fastq_2 ] ]
    gff       = BAKTA_BAKTA.out.gff
    scaffolds = FILTER_ASSEMBLY.out.filtered

    versions = ch_versions                     // channel: [ versions.yml ]
}

