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
    sample_data
        .map{ [[id: it.sample_id, single_end: it.single_end], it.paths, it.sequence_type, it.kingdom] }
        .unique()
        .branch { meta, paths, type, kingdom ->
            short_prokaryote:    (type == "illumina" || type == "bgiseq") && kingdom == "bacteria"
                return [meta, paths]
            nanopore_prokaryote: type == "nanopore" && kingdom == "bacteria"
                return [meta, paths]
            pacbio_prokaryote:   type == "pacbio" && kingdom == "bacteria"
                return [meta, paths]
            short_eukaryote:     (type == "illumina" || type == "bgiseq") && kingdom != "bacteria"
                return [meta, paths]
            nanopore_eukaryote:  type == "nanopore" && kingdom != "bacteria"
                return [meta, paths]
            pacbio_eukaryote:    type == "pacbio" && kingdom != "bacteria"
                return [meta, paths]
            other:               true
                return [meta, paths]
        }
        .set { filtered_input }

    FASTP( filtered_input.short_prokaryote, [], false, false )
    versions = versions.mix(FASTP.out.versions)

    SPADES(
        FASTP.out.reads.map{ sample_meta, read_paths -> [sample_meta, read_paths, [], []] },
        [], // val yml
        []  // val hmm
    )
    versions = versions.mix(SPADES.out.versions)

    FLYE_NANOPORE (
        filtered_input.nanopore_prokaryote.mix(filtered_input.nanopore_eukaryote),
        "--nano-hq"
    )

    FLYE_PACBIO (
        filtered_input.pacbio_prokaryote.mix(filtered_input.pacbio_eukaryote),
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

    // Warn if a sample was not assembled
    not_assembled_warnings = sample_data
        .map { [[id: it.sample_id], it] }
        .combine(filtered_input.other, by: 0)
        .map{ sample_meta, sample_data, paths ->
            [sample_meta, [id: sample_data.report_group_ids], null, "GENOME_ASSEMBLY", "WARNING", "Sample not assembled because no assemblier was configured to handle this combination of taxon and sequencing technology"]
        }
    messages = messages.mix(not_assembled_warnings)

    emit:
    reads     = FASTP.out.reads           // channel: [ val(meta), [reads] ]
    scaffolds = filtered_assembly
    quast     = QUAST.out.results
    versions  = versions
    messages  = messages    // meta, group_meta, ref_meta, workflow, level, message
}

