include { FASTP                 } from '../../../modules/nf-core/fastp'
include { SPADES                } from '../../../modules/nf-core/spades'
include { QUAST                 } from '../../../modules/nf-core/quast'
include { FLYE as FLYE_NANOPORE } from '../../../modules/nf-core/flye'
include { FLYE as FLYE_PACBIO   } from '../../../modules/nf-core/flye'

workflow GENOME_ASSEMBLY {

    take:
    sample_data

    main:

    versions = channel.empty()
    messages = channel.empty()
    parsed_sample_data = sample_data
        .map{ [[id: it.sample_id, single_end: it.single_end, domain: it.domain, type: it.sequence_type], [id: it.report_group_ids], it.paths] }
    parsed_sample_data
        .map{ sample_meta, report_meta, read_paths -> [sample_meta, read_paths]}
        .unique()
        .branch { meta, paths->
            short_prokaryote:    (meta.type == "illumina" || meta.type == "bgiseq") && meta.domain == "Bacteria"
            nanopore_prokaryote: meta.type == "nanopore" && meta.domain == "Bacteria"
            pacbio_prokaryote:   meta.type == "pacbio" && meta.domain == "Bacteria"
            short_eukaryote:     (meta.type == "illumina" || meta.type == "bgiseq") && meta.domain != "Bacteria"
            nanopore_eukaryote:  meta.type == "nanopore" && meta.domain != "Bacteria"
            pacbio_eukaryote:    meta.type == "pacbio" && meta.domain != "Bacteria"
            other:               true
        }
        .set { filtered_input }

    spades_input = filtered_input.short_prokaryote
        .mix(filtered_input.short_eukaryote)
    fastp_input = spades_input
        .map{ sample_meta, read_paths ->    // If there are both single and paired in reads, just use the paired end reads
            [sample_meta, read_paths.size() <= 2 ? read_paths : read_paths.findAll { it ==~ /.+_[12]\..+$/ }, [] ]
        }
    FASTP( fastp_input, false, false, false )

    // Check for samples with too few reads after quality control
    filtered_reads = FASTP.out.json
        .splitJson(path: 'summary.after_filtering')
        .filter{ sample_meta, json ->
            json.key == 'total_bases'
        }
        .combine(FASTP.out.reads, by: 0)
        .branch{ sample_meta, json, read_paths ->
            pass: json.value.toBigInteger() >= params.min_bases_to_assemble.toBigInteger()
                return [sample_meta, read_paths]

            fail: true
                return [sample_meta, read_paths, json.value]
        }

    filtered_reads_warnings = filtered_reads.fail
        .combine(parsed_sample_data, by: 0)
        .map { sample_meta, read_paths1, base_count, report_meta, read_paths2 ->
            [sample_meta, report_meta, null, "GENOME_ASSEMBLY", "WARNING", "After quality filtering, sample reads consist of ${base_count} bases, which is less than the minimum of ${params.min_bases_to_assemble} defined by the option `min_bases_to_assemble` and therefore will not be assembled."]
        }
    messages = messages.mix(filtered_reads_warnings)

    SPADES(
        filtered_reads.pass.map{ sample_meta, read_paths -> [sample_meta, read_paths, [], []] },
        [], // val yml
        []  // val hmm
    )
    versions = versions.mix(SPADES.out.versions)

    // Warn about any failed Spades assemblies
    spades_warnings = spades_input
        .join(SPADES.out.scaffolds, remainder: true)
        .filter { sample_meta, read_paths, scaffolds ->
            ! scaffolds
        }
        .combine(parsed_sample_data, by: 0)
        .map { sample_meta, read_paths, scaffolds, report_meta, read_paths2 ->
            [sample_meta, report_meta, null, "GENOME_ASSEMBLY", "WARNING", "Sample could not be assembled, possibly due to short read lengh or low quality. Check Spades' logs for more details."]
        }
    messages = messages.mix(spades_warnings)

    FLYE_NANOPORE (
        filtered_input.nanopore_prokaryote.mix(filtered_input.nanopore_eukaryote),
        "--nano-raw"
    )

    FLYE_PACBIO (
        filtered_input.pacbio_prokaryote.mix(filtered_input.pacbio_eukaryote),
        "--pacbio-raw"
    )

    assemblies = SPADES.out.scaffolds
        .mix(FLYE_NANOPORE.out.fasta)
        .mix(FLYE_PACBIO.out.fasta)
        .map { sample_meta, path ->  // remove the "single_end" in the sample meta data so that it is just the ID like most of the pipeline
            [[id: sample_meta.id], path]
        }
        .unique()
    QUAST (
        assemblies, [[], []], [[],[]]
    )

    // Warn if a sample was not assembled
    not_assembled_warnings = sample_data
        .map { [[id: it.sample_id], it] }
        .combine(filtered_input.other, by: 0)
        .map{ sample_meta, sample_data, paths ->
            [sample_meta, [id: sample_data.report_group_ids], null, "GENOME_ASSEMBLY", "WARNING", "Sample not assembled because no assemblier was configured to handle this combination of taxon and sequencing technology"]
        }
    messages = messages.mix(not_assembled_warnings)

    emit:
    reads      = FASTP.out.reads.map{sample_meta, reads -> [[id: sample_meta.id], reads]} // strip off extra sample metadata to make joins easier downstream
    fastp_json = FASTP.out.json.map{sample_meta, json -> [[id: sample_meta.id], json]} // strip off extra sample metadata to make joins easier downstream
    scaffolds  = assemblies
    quast      = QUAST.out.results
    versions   = versions
    messages   = messages    // meta, group_meta, ref_meta, workflow, level, message
}
