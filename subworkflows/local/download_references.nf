include { FIND_ASSEMBLIES                           } from '../../modules/local/find_assemblies'
include { MERGE_ASSEMBLIES                          } from '../../modules/local/merge_assemblies'
include { PICK_ASSEMBLIES                           } from '../../modules/local/pick_assemblies'
include { DOWNLOAD_ASSEMBLIES                       } from '../../modules/local/download_assemblies'
include { MAKE_GFF_WITH_FASTA                       } from '../../modules/local/make_gff_with_fasta'
include { SOURMASH_SKETCH as SOURMASH_SKETCH_GENOME } from '../../modules/nf-core/sourmash/sketch/main'
include { KHMER_TRIMLOWABUND                        } from '../../modules/local/khmer_trimlowabund'

workflow DOWNLOAD_REFERENCES {

    take:
    ch_species  // channel: [val(meta), file(taxa)]
    ch_genera  // channel: [val(meta), file(taxa)]
    ch_families  // channel: [val(meta), file(taxa)]
    ch_input // channel: [val(meta), [file(fastq)], val(sra), val(ref_meta), file(reference), val(ref_acc), val(group_meta)]

    main:
    ch_versions = Channel.empty()

    ch_all_families = ch_families
        .map {it[1]}
        .splitText()
        .map { it.replace('\n', '') }
        .collect()
        .toSortedList()
        .flatten()
        .unique()

    FIND_ASSEMBLIES ( ch_all_families )
    ch_versions = ch_versions.mix(FIND_ASSEMBLIES.out.versions.toSortedList().map{it[0]})

    PICK_ASSEMBLIES (
        ch_families
            .join(ch_genera)
            .join(ch_species),
        FIND_ASSEMBLIES.out.stats
            .map { it[1] }
            .toSortedList()
    )

    // Make channel with all unique assembly IDs
    user_acc_list = ch_input
        .map { [it[3], it[5]] } // ref_meta, ref_acc
        .distinct()
    ch_assembly_ids = PICK_ASSEMBLIES.out.id_list
        .map {it[1]}
        .splitText()
        .map { it.replace('\n', '') }
        .filter { it != '' }
        .map { it.split('\t') }
        .map { [[id: it[0]], it[1]] }
        .unique()
        .join(user_acc_list, by:1, remainder: true) // this is used to provide something for the following filter to work
        .filter {it[2] == null} // remove any user-defined accession numbers that have already been downloaded
        .map { it[1..0] }
        .unique()
    DOWNLOAD_ASSEMBLIES ( ch_assembly_ids )
    ch_versions = ch_versions.mix(DOWNLOAD_ASSEMBLIES.out.versions.toSortedList().map{it[0]})

    // Add sequence to the gff
    MAKE_GFF_WITH_FASTA (
        DOWNLOAD_ASSEMBLIES.out.sequence
            .join(DOWNLOAD_ASSEMBLIES.out.gff)
    )

    SOURMASH_SKETCH_GENOME ( DOWNLOAD_ASSEMBLIES.out.sequence )
    ch_versions = ch_versions.mix(SOURMASH_SKETCH_GENOME.out.versions.toSortedList().map{it[0]})

    genome_ids = PICK_ASSEMBLIES.out.id_list
        .splitText(elem: 1)
        .map { [[id: it[1].replace('\n', '').split('\t')[0]], it[0]] } // [ val(ref_meta), val(meta) ]

    emit:
    assem_samp_combos = genome_ids                        // [ val(ref_meta), val(meta) ] for each assembly/sample combination
    sequence   = DOWNLOAD_ASSEMBLIES.out.sequence         // [ val(ref_meta), file(fna) ] for each assembly
    gff        = MAKE_GFF_WITH_FASTA.out.gff              // [ val(ref_meta), file(gff) ] for each assembly
    signatures = SOURMASH_SKETCH_GENOME.out.signatures    // [ val(ref_meta), file(signature) ] for each assembly
    stats      = PICK_ASSEMBLIES.out.merged_stats.first() // [ file(stats) ]
    versions   = ch_versions                              // [ versions.yml ]
}
