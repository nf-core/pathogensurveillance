include { FASTP                 } from '../../../modules/nf-core/fastp'
include { SPADES                } from '../../../modules/nf-core/spades'
include { QUAST                 } from '../../../modules/nf-core/quast'
include { FLYE as FLYE_NANOPORE } from '../../../modules/nf-core/flye'
include { FLYE as FLYE_PACBIO   } from '../../../modules/nf-core/flye'

def readsSizeGB(read_paths) {
    if (read_paths instanceof List) {
        return read_paths.collect { path -> path.size() }.sum() / (1024 * 1024 * 1024)
    } else {
        return read_paths.size / (1024 * 1024 * 1024)
    }
}

workflow GENOME_ASSEMBLY {

    take:
    sample_data

    main:

    versions = channel.empty()
    messages = channel.empty()
    parsed_sample_data = sample_data
        .map{ sample_meta -> [[id: sample_meta.sample_id, single_end: sample_meta.single_end, domain: sample_meta.domain, type: sample_meta.sequence_type], [id: sample_meta.report_group_ids], sample_meta.paths] }
    all_samples = parsed_sample_data
        .map{ sample_meta, report_meta, read_paths ->
            [sample_meta + [reads_size_gb: readsSizeGB(read_paths)], read_paths]
        }
        .unique()

    short_reads = all_samples.filter { meta, paths -> meta.type == "illumina" || meta.type == "bgiseq" }
    nanopore_reads = all_samples.filter { meta, paths -> meta.type == "nanopore" }
    pacbio_reads = all_samples.filter { meta, paths -> meta.type == "pacbio" }
    other_reads = all_samples.filter { meta, paths ->
        meta.type != "illumina" && meta.type != "bgiseq" && meta.type != "nanopore" && meta.type != "pacbio"
    }

    spades_input = short_reads
    fastp_input = spades_input
        .map{ sample_meta, read_paths ->    // If there are both single and paired in reads, just use the paired end reads
            [sample_meta, read_paths.size() <= 2 ? read_paths : read_paths.findAll { path -> path ==~ /.+_[12]\..+$/ }, [] ]
        }
    FASTP( fastp_input, false, false, false )

    // Check for samples with too few reads after quality control
    fastp_combined = FASTP.out.json
        .splitJson(path: 'summary.after_filtering')
        .filter{ sample_meta, json ->
            json.key == 'total_bases'
        }
        .combine(FASTP.out.reads, by: 0)

    filtered_reads_pass = fastp_combined
        .filter { sample_meta, json, read_paths -> json.value.toBigInteger() >= params.min_bases_to_assemble.toBigInteger() }
        .map { sample_meta, json, read_paths -> [sample_meta, read_paths] }

    filtered_reads_fail = fastp_combined
        .filter { sample_meta, json, read_paths -> json.value.toBigInteger() < params.min_bases_to_assemble.toBigInteger() }
        .map { sample_meta, json, read_paths -> [sample_meta, read_paths, json.value] }

    filtered_reads_warnings = filtered_reads_fail
        .combine(parsed_sample_data, by: 0)
        .map { sample_meta, read_paths1, base_count, report_meta, read_paths2 ->
            [sample_meta, report_meta, null, "GENOME_ASSEMBLY", "WARNING", "After quality filtering, sample reads consist of ${base_count} bases, which is less than the minimum of ${params.min_bases_to_assemble} defined by the option `min_bases_to_assemble` and therefore will not be assembled."]
        }
    messages = messages.mix(filtered_reads_warnings)

    SPADES(
        filtered_reads_pass.map{ sample_meta, read_paths -> [sample_meta, read_paths, [], []] },
        [], // val yml
        []  // val hmm
    )

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

    FLYE_NANOPORE (nanopore_reads, "--nano-raw")
    FLYE_PACBIO (pacbio_reads, "--pacbio-raw")

    // Warn about any failed Flye assemblies
    flye_warnings = FLYE_NANOPORE.out.fasta
        .mix(FLYE_PACBIO.out.fasta)
        .filter { sample_meta, fasta ->
            ! fasta || fasta.size() == 0
        }
        .combine(parsed_sample_data, by: 0)
        .map { sample_meta, fasta, report_meta, read_paths2 ->
            [sample_meta, report_meta, null, "GENOME_ASSEMBLY", "WARNING", "Sample could not be assembled. Check Flye logs for more details."]
        }
    messages = messages.mix(flye_warnings)

    assemblies = SPADES.out.scaffolds
        .mix(FLYE_NANOPORE.out.fasta)
        .mix(FLYE_PACBIO.out.fasta)
        .filter { sample_meta, path -> path && path.size() > 0 }  // Filter out empty assemblies
        .map { sample_meta, path ->  // remove the "single_end" in the sample meta data so that it is just the ID like most of the pipeline
            [[id: sample_meta.id], path]
        }
        .unique()
    QUAST (
        assemblies, [[], []], [[],[]]
    )

    // Warn if a sample was not assembled
    not_assembled_warnings = sample_data
        .map { sample_meta -> [[id: sample_meta.sample_id], sample_meta] }
        .combine(other_reads, by: 0)
        .map{ sample_meta, my_sample_data, paths ->
            [sample_meta, [id: my_sample_data.report_group_ids], null, "GENOME_ASSEMBLY", "WARNING", "Sample not assembled because no assemblier was configured to handle this combination of taxon and sequencing technology"]
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
