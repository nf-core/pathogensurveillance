include { FIND_ASSEMBLIES       } from '../../modules/local/findassemblies'
include { MERGE_ASSEMBLIES      } from '../../modules/local/mergeassemblies'
include { PICK_ASSEMBLIES       } from '../../modules/local/pickassemblies'
include { DOWNLOAD_ASSEMBLIES   } from '../../modules/local/downloadassemblies'
include { MAKEGFFWITHFASTA      } from '../../modules/local/makegffwithfasta'

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
        ch_families.join(ch_genera).join(ch_species),
        MERGE_ASSEMBLIES.out.merged_stats
    )

    // Make channel with all unique assembly IDs 
    ch_assembly_ids = PICK_ASSEMBLIES.out.id_list
        .map {it[1]}                                                            
        .splitText()                                                            
        .map { it.replace('\n', '') }                                           
        .collect()                                                              
        .toSortedList()                                                         
        .flatten()                                                              
        .unique()
    DOWNLOAD_ASSEMBLIES ( ch_assembly_ids )
    
    MAKEGFFWITHFASTA ( DOWNLOAD_ASSEMBLIES.out.sequence.join(DOWNLOAD_ASSEMBLIES.out.gff) )

    ch_samp_ids = PICK_ASSEMBLIES.out.id_list
        .splitText(elem: 1)
        .map { [it[1].replace('\n', ''), it[0]] }

    ch_samp_seq = DOWNLOAD_ASSEMBLIES.out.sequence
        .cross(ch_samp_ids)           // The cross/map combo does a full join
        .map { [it[1][1], it[0][1]] } //   (join operator is an inner join).
        .groupTuple()

    ch_samp_gff = DOWNLOAD_ASSEMBLIES.out.gff                            
        .cross(ch_samp_ids)           // The cross/map combo does a full join   
        .map { [it[1][1], it[0][1]] } //   (join operator is an inner join).    
        .groupTuple().view()                                                         

    emit:
    stats           = FIND_ASSEMBLIES.out.stats        // channel: [ val(taxon), file(stats) ]
    versions        = ch_versions                      // channel: [ versions.yml ]
}
