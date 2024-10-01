include { FASTP                 } from '../../modules/nf-core/fastp/main'
include { SPADES                } from '../../modules/nf-core/spades/main'
include { FILTER_ASSEMBLY       } from '../../modules/local/filter_assembly'
include { QUAST                 } from '../../modules/local/quast.nf'
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
        .map{ [[id: it.sample_id, single_end: it.single_end], it.paths, it.sequence_type] }
        .unique()

    shortreads = filtered_input
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

    nanopore = filtered_input
        .filter{ sample_meta, read_paths, seq_type -> seq_type == "nanopore"}
        .map{ sample_meta, read_paths, seq_type -> [sample_meta, read_paths] }
    FLYE_NANOPORE (
        nanopore,
        "--nano-hq"
    )

    pacbio = filtered_input
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
        .map { sample_meta, path ->  // remove the "single_end" in the sample meta data so that it is just the ID like most of the pipeline
            [[id: sample_meta.id], path]
        }
        .unique()
    QUAST (
        filtered_assembly
            .map { sample_meta, assembly ->
                [sample_meta, assembly, [], []]
            }
    )
    versions = versions.mix(QUAST.out.versions)

    emit:
    reads     = FASTP.out.reads           // channel: [ val(meta), [reads] ]
    scaffolds = filtered_assembly
    quast     = QUAST.out.results
    versions  = versions
    messages  = messages    // meta, group_meta, ref_meta, workflow, level, message
}

