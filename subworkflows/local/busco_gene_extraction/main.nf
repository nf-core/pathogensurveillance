include { BUSCO_BUSCO        } from '../../../modules/nf-core/busco/busco'
include { BUSCO_DOWNLOAD     } from '../../../modules/nf-core/busco/download'
include { FORMAT_BUSCO_GENES } from '../../../modules/local/format_busco_genes'

workflow BUSCO_GENE_EXTRACTION {

    take:
    sample_data       // already domain-filtered and restricted to samples with assemblies
    sample_assemblies // sample_id, assembly_path
    selected_ref_data // ref_meta, report_meta, ref_path, ref_name, gff_path

    main:

    messages = channel.empty()

    // Combine selected reference data with sample metadata and assembled genomes
    selected_ref_input = selected_ref_data
        .map { ref_meta, report_meta, ref_path, ref_name, gff_path ->
            [ref_meta, report_meta, ref_path]
        }

    busco_input = sample_data
        .map{ sample_meta -> [[id: sample_meta.sample_id], [id: sample_meta.report_group_ids]] }
        .combine(sample_assemblies, by: 0)
        .mix(selected_ref_input)
        .unique()

    // Download BUSCO datasets
    def busco_cache = file("${params.data_dir}/busco_db/busco_downloads")
    if (busco_cache.exists()) {
        busco_download = Channel.fromPath("${params.data_dir}/busco_db/busco_downloads").first()
    } else {
        BUSCO_DOWNLOAD ( channel.from( "eukaryota_odb10" ) )
        busco_download = BUSCO_DOWNLOAD.out.download_dir.first() // .first() is needed to convert the queue channel to a value channel so it can be used multiple times.
    }

    // Extract BUSCO genes for all unique reference genomes used in any sample/group
    BUSCO_BUSCO (
        busco_input
            .map{ meta, report_meta, path ->
                [meta, path]
            }
            .unique(),
        "genome",
        "eukaryota_odb10",
        busco_download,
        [],
        true
    )

    // Combine BUSCO output by report group and convert to unified gene_sequences format
    sorted_busco_data = busco_input
        .combine(BUSCO_BUSCO.out.busco_dir, by: 0)
        .map{ meta, report_meta, path, busco_dir ->
            [report_meta, busco_dir]
        }
        .groupTuple(by: 0, sort: 'hash')

    FORMAT_BUSCO_GENES(sorted_busco_data)

    emit:
    gene_sequences = FORMAT_BUSCO_GENES.out.gene_sequences
        .map { group_meta, gene_seqs -> [group_meta, gene_seqs, 'busco', 'BUSCO_GENE_EXTRACTION'] }
    messages       = messages

}
