include { BUSCO          } from '../../modules/nf-core/busco/main'
include { BUSCO_DOWNLOAD } from '../../modules/local/busco_download'
include { READ2TREE      } from '../../modules/local/read2tree/main'

workflow BUSCO_PHYLOGENY {

    take:
    input // val(meta), [file(fastq)], val(group_meta), [val(ref_meta)], val(kingdom), val(depth)
    ref_seqs // val(ref_meta), file(fna)

    main:

    ch_versions = Channel.empty()
    messages = Channel.empty()

    // Only process Eukaryotic samples
    input_filtered = input
        .filter { it[4] == "Eukaryota" }

    // Get all reference genomes in any group/sample
    unique_ref_fnas = input_filtered
        .map { [it[3]] }
        .flatten()
        .unique()
        .join(ref_seqs)  // val(ref_meta), file(fna)

    // Download BUSCO datasets
    BUSCO_DOWNLOAD ( Channel.from( "eukaryota_odb10" ) )

    // Extract BUSCO genes for all unique reference genomes used in any sample/group
    BUSCO (
        unique_ref_fnas,
        "genome",
        "eukaryota_odb10",
        BUSCO_DOWNLOAD.out.download_dir.first(), // .first() is needed to convert the queue channel to a value channel so it can be used multiple times.
        [] )

    // Create Read2tree database

    // Run Read2tree
    // READ2TREE ( )


    emit:
    versions   = ch_versions // versions.yml
    messages   = messages    // meta, group_meta, ref_meta, workflow, level, message

}
