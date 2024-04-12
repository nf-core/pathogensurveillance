include { BUSCO          } from '../../modules/nf-core/busco/main'
include { BUSCO_DOWNLOAD } from '../../modules/local/busco_download'
include { READ2TREE      } from '../../modules/local/read2tree/main'
include { R2TF           } from '../../modules/local/r2tf'
include { R2TDIR         } from '../../modules/local/r2tdir'
include { R2TBIN         } from '../../modules/local/r2tbin'

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
    R2TF (
        BUSCO.out.busco_dir
	)

    // Create directories for Read2Tree
    R2TDIR (
        R2TF.out.output
            .map {it[1]}
            .collect()
        )
    // Bin busco files by genes
    R2TBIN (
        R2TDIR.out.markers
    )
    // Run Read2tree
    // READ2TREE ( )


    emit:
    versions   = ch_versions // versions.yml
    messages   = messages    // meta, group_meta, ref_meta, workflow, level, message

}
