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
    sample_data

    main:

    versions = Channel.empty()
    messages = Channel.empty()
    filtered_input = sample_data
        .filter {it.kingdom == "Bacteria"}
        .map{ [[id: it.sample_id], it.paths, it.kingdom, it.sendsketch_depth, it.sequence_type] }
        .unique()
    reads = filtered_input
        .map { sample_meta, read_paths, kingdom, depth, seq_type ->
            [sample_meta, read_paths, depth]
        }
        .unique()

    // Subset sample reads to increase speed of following steps
    SUBSET_READS (
        reads,
        params.sketch_max_depth
    )
    versions = versions.mix(SUBSET_READS.out.versions)
    subset_reads = SUBSET_READS.out.reads
        .join(filtered_input)
        .map { sample_meta, subset_read_paths, read_paths, kingdom, depth, seq_type ->
            [sample_meta, subset_read_paths, seq_type]
        }

    shortreads = subset_reads
        .filter{ sample_meta, read_paths, seq_type -> seq_type == "illumina" || seq_type == "bgiseq" }
        .map{ sample_meta, read_paths, seq_type -> [sample_meta, read_paths] }
    FASTP( shortreads, [], false, false )
    versions = versions.mix(FASTP.out.versions)

    SPADES(
        FASTP.out.reads.map{ sample_meta, read_paths -> [sample_meta, read_paths, [], []] },
        [], // val yml
        []  // val hmm
    )
    versions = versions.mix(SPADES.out.versions)

    nanopore = subset_reads
        .filter{ sample_meta, read_paths, seq_type -> seq_type == "nanopore"}
        .map{ sample_meta, read_paths, seq_type -> [sample_meta, read_paths] }
    FLYE_NANOPORE (
        nanopore,
        "--nano-hq"
    )

    pacbio = subset_reads
        .filter{ sample_meta, read_paths, seq_type -> seq_type == "pacbio"}
        .map{ sample_meta, read_paths, seq_type -> [sample_meta, read_paths] }
    FLYE_PACBIO (
        pacbio,
        "--pacbio-hifi"
    )

    FILTER_ASSEMBLY (
        SPADES.out.scaffolds
    )
    versions = versions.mix(FILTER_ASSEMBLY.out.versions)

    filtered_assembly = FILTER_ASSEMBLY.out.filtered
        .mix(FLYE_NANOPORE.out.fasta)
        .mix(FLYE_PACBIO.out.fasta)
        .unique()
    QUAST (
        filtered_assembly
            .map { sample_meta, assembly ->
                [sample_meta, assembly, [], []]
            }
    )
    versions = versions.mix(QUAST.out.versions)

    // Download the bakta database if needed
    //   Based on code from the bacass nf-core pipeline using the MIT license: https://github.com/nf-core/bacass
    if (params.bakta_db) {
        if (params.bakta_db.endsWith('.tar.gz')) {
            bakta_db_tar = Channel.fromPath(params.bakta_db).map{ [ [id: 'baktadb'], it] }
            UNTAR( bakta_db_tar )
            bakta_db = UNTAR.out.untar.map{ meta, db -> db }.first()
            versions = versions.mix(UNTAR.out.versions)
        } else {
            bakta_db = Channel.fromPath(params.bakta_db).first()
        }
    } else if (params.download_bakta_db) {
        BAKTA_BAKTADBDOWNLOAD()
        bakta_db = BAKTA_BAKTADBDOWNLOAD.out.db
        versions = versions.mix(BAKTA_BAKTADBDOWNLOAD.out.versions)
    }

    // Run bakta
    BAKTA_BAKTA (
        filtered_assembly, // Genome assembly
        bakta_db, // Bakta database
        [], // proteins (optional)
        [] // prodigal_tf (optional)
    )
    versions = versions.mix(BAKTA_BAKTA.out.versions)

    emit:
    reads     = FASTP.out.reads           // channel: [ val(meta), [reads] ]
    gff       = BAKTA_BAKTA.out.gff
    scaffolds = filtered_assembly
    quast     = QUAST.out.results
    versions  = versions
    messages  = messages    // meta, group_meta, ref_meta, workflow, level, message
}

