include { FASTP                 } from '../../modules/nf-core/fastp/main'
include { SPADES                } from '../../modules/nf-core/spades/main'
include { FILTER_ASSEMBLY       } from '../../modules/local/filter_assembly'
include { QUAST                 } from '../../modules/local/quast.nf'
include { BAKTA_BAKTA           } from '../../modules/nf-core/bakta/bakta/main'
include { BAKTA_BAKTADBDOWNLOAD } from '../../modules/nf-core/bakta/baktadbdownload/main'
include { SUBSET_READS          } from '../../modules/local/subset_reads'
include { UNTAR                 } from '../../modules/nf-core/untar/main'
include { FLYE as FLYE_NANOPORE } from '../../modules/nf-core/flye/main'
include { FLYE as FLYE_PACBIO   } from '../../modules/nf-core/flye/main'

workflow GENOME_ASSEMBLY {

    take:
    ch_input // meta, [reads], ref_meta, reference, group_meta, kingdom, depth

    main:

    ch_versions = Channel.empty()
    ch_input_filtered = ch_input
        .filter { it[5] == "Bacteria" }
    ch_reads = ch_input_filtered
        .map { [it[0], it[1], it[8]] }
        .unique()

    // Subset sample reads to increase speed of following steps
    SUBSET_READS (
        ch_reads,
        params.sketch_max_depth
    )
    ch_versions = ch_versions.mix(SUBSET_READS.out.versions.first())
    subset_reads = SUBSET_READS.out.reads
        .join(ch_input_filtered) // meta, [subset_reads], [reads], ref_meta, reference, group_meta, kingdom, depth

    shortreads = subset_reads
        .filter {it[0].id == "illumina"}
        .map { it[0..1] } // meta, [subset_reads]
    FASTP (shortreads, [], false, false )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    SPADES (
        FASTP.out.reads.map { [it[0], it[1], [], []] }, // val(meta), path(illumina), path(pacbio), path(nanopore)
        [], // val yml
        []  // val hmm
    )
    ch_versions = ch_versions.mix(SPADES.out.versions.first())

    nanopore = subset_reads
        .filter {it[0].reads_type == "nanopore"}
        .map { it[0..1] }
    FLYE_NANOPORE (
        nanopore,
        "--nano-hq"
    )

    pacbio = subset_reads
        .filter {it[0].reads_type == "pacbio"}
        .map { it[0..1] }
    FLYE_PACBIO (
        pacbio,
        "--pacbio-hifi"
    )

    FILTER_ASSEMBLY (
        SPADES.out.scaffolds
    )
    ch_versions = ch_versions.mix(FILTER_ASSEMBLY.out.versions.first())

    filtered_assembly = FILTER_ASSEMBLY.out.filtered
        .mix(FLYE_NANOPORE.out.fasta)
        .mix(FLYE_PACBIO.out.fasta)
    ch_ref_grouped = ch_input_filtered
        .combine(filtered_assembly, by: 0) // meta, [reads], ref_meta, reference, group_meta, kingdom, depth, filt_assemb
        .groupTuple(by: 2) // meta, [reads], ref_meta, reference, group_meta, kingdom, depth, filt_assemb
        .map { [it[2], it[7].sort().unique(), it[3].sort()[0] ?: [], []] } // ref_meta, assembly, reference, gff
    QUAST ( ch_ref_grouped )
    ch_versions = ch_versions.mix(QUAST.out.versions.first())

    // Download the bakta database if needed
    //   Based on code from the bacass nf-core pipeline using the MIT license: https://github.com/nf-core/bacass
    if (params.bakta_db) {
        if (params.bakta_db.endsWith('.tar.gz')) {
            bakta_db_tar = Channel.fromPath(params.bakta_db).map{ [ [id: 'baktadb'], it] }
            UNTAR( bakta_db_tar )
            bakta_db = UNTAR.out.untar.map{ meta, db -> db }.first()
            ch_versions = ch_versions.mix(UNTAR.out.versions)
        } else {
            bakta_db = Channel.fromPath(params.bakta_db).first()
        }
    } else if (params.download_bakta_db){
        BAKTA_BAKTADBDOWNLOAD()
        bakta_db  = BAKTA_BAKTADBDOWNLOAD.out.db
        ch_versions = ch_versions.mix(BAKTA_BAKTADBDOWNLOAD.out.versions)
    }

    // Run bakta
    BAKTA_BAKTA (
        filtered_assembly, // Genome assembly
        bakta_db, // Bakta database
        [], // proteins (optional)
        [] // prodigal_tf (optional)
    )
    ch_versions = ch_versions.mix(BAKTA_BAKTA.out.versions.first())

    emit:
    reads     = FASTP.out.reads           // channel: [ val(meta), [reads] ]
    gff       = BAKTA_BAKTA.out.gff
    scaffolds = filtered_assembly
    quast     = QUAST.out.results

    versions = ch_versions                     // channel: [ versions.yml ]
}

