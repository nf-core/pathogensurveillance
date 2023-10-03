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

    MERGE_ASSEMBLIES (
       FIND_ASSEMBLIES.out.stats
           .map { it[1] }
           .collect()
    )
    
    PICK_ASSEMBLIES (
        ch_families
            .join(ch_genera)
            .join(ch_species),
        MERGE_ASSEMBLIES.out.merged_stats
    )

    // Make channel with all unique assembly IDs 
    ch_assembly_ids = PICK_ASSEMBLIES.out.id_list
        .map {it[1]}                                                            
        .splitText()                                                            
        .map { it.replace('\n', '') }
        .filter { it != '' }
        .collect()                                                              
        .toSortedList()                                                         
        .flatten()                                                              
        .unique()
    DOWNLOAD_ASSEMBLIES ( ch_assembly_ids )
    ch_versions = ch_versions.mix(DOWNLOAD_ASSEMBLIES.out.versions.toSortedList().map{it[0]})
    
    MAKE_GFF_WITH_FASTA ( DOWNLOAD_ASSEMBLIES.out.sequence.join(DOWNLOAD_ASSEMBLIES.out.gff) )
   
    SOURMASH_SKETCH_GENOME (
        DOWNLOAD_ASSEMBLIES.out.sequence
        .map { [[id: it[0]], it[1]] }
    )
    ch_versions = ch_versions.mix(SOURMASH_SKETCH_GENOME.out.versions.toSortedList().map{it[0]})

    genome_ids = PICK_ASSEMBLIES.out.id_list
        .splitText(elem: 1)
        .map { [it[1].replace('\n', ''), it[0]] } // [ val(genome_id), val(meta) ]


    emit:
    assem_samp_combos = genome_ids                            // [ val(genome_id), val(meta) ] for each assembly/sample combination
    sequence   = DOWNLOAD_ASSEMBLIES.out.sequence      // [ val(genome_id), file(fna) ] for each assembly
    gff        = MAKE_GFF_WITH_FASTA.out.gff              // [ val(genome_id), file(gff) ] for each assembly
    signatures = SOURMASH_SKETCH_GENOME.out.signatures // [ val(genome_id), file(signature) ] for each assembly
    stats      = MERGE_ASSEMBLIES.out.merged_stats     // [ file(stats) ]
    versions   = ch_versions                           // [ versions.yml ]
}
