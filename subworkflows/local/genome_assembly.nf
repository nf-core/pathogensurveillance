include { FASTP           } from '../../modules/nf-core/fastp/main'
include { SPADES          } from '../../modules/nf-core/spades/main'
include { FILTER_ASSEMBLY } from '../../modules/local/filter_assembly'
include { QUAST           } from '../../modules/local/quast.nf'
include { BAKTA_BAKTA     } from '../../modules/nf-core/bakta/bakta/main'
include { SUBSET_READS    } from '../../modules/local/subset_reads'                 

workflow GENOME_ASSEMBLY {

    take:
    ch_input // channel: [ val(meta), [fastq_1, fastq_2], val(ref_meta), file(reference), val(group_meta), val(kingdom), val(depth) ]

    main:

    ch_versions = Channel.empty()
    ch_input_filtered = ch_input
        .filter { it[5] == "Bacteria" }
    ch_reads = ch_input_filtered
        .map { it[0..1] + [it[6]] }
        .unique()
    
    // Subset sample reads to increase speed of following steps                 
    SUBSET_READS (
        ch_reads,                                                              
        params.sketch_max_depth                                                 
    )                                                                           
    ch_versions = ch_versions.mix(SUBSET_READS.out.versions.first())

    FASTP ( SUBSET_READS.out.reads, [], false, false )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    SPADES (
        FASTP.out.reads.map { [it[0], it[1], [], []] }, // val(meta), path(illumina), path(pacbio), path(nanopore)
        [], // val yml
        []  // val hmm
    )
    ch_versions = ch_versions.mix(SPADES.out.versions.first())

    FILTER_ASSEMBLY (
        SPADES.out.scaffolds
    )
    ch_versions = ch_versions.mix(FILTER_ASSEMBLY.out.versions.first())

    ch_ref_grouped = ch_input_filtered
        .combine(FILTER_ASSEMBLY.out.filtered, by: 0)
        .groupTuple(by: 2) // [val(meta)], [[fastq_1, fastq_2]], val(ref_meta), [file(reference)], [val(group_meta)], [val(kingdom)], val(depth), file(assembly)]
        .map { [it[2], it[7].sort(), it[3].sort()[0], []] } // ref_meta, assembly, reference
    QUAST (
        ch_ref_grouped,
        true, // use_fasta
        false // use_gff
    )
    ch_versions = ch_versions.mix(QUAST.out.versions.first())

    ch_bakta_db = Channel.value("$projectDir/assets/bakta_db/db-light")
    BAKTA_BAKTA (
        FILTER_ASSEMBLY.out.filtered, // Genome assembly
        ch_bakta_db, // Bakta database
        [], // proteins (optional)
        [] // prodigal_tf (optional)
    )
    ch_versions = ch_versions.mix(BAKTA_BAKTA.out.versions.first())

    emit:
    reads     = FASTP.out.reads           // channel: [ val(meta), [ fastq_1, fastq_2 ] ]
    gff       = BAKTA_BAKTA.out.gff
    scaffolds = FILTER_ASSEMBLY.out.filtered
    quast     = QUAST.out.results

    versions = ch_versions                     // channel: [ versions.yml ]
}

