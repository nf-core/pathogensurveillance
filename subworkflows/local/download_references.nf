include { FIND_ASSEMBLIES       } from '../../modules/local/findassemblies'
include { MERGE_ASSEMBLIES      } from '../../modules/local/mergeassemblies'
include { PICK_ASSEMBLIES       } from '../../modules/local/pickassemblies'
include { DOWNLOAD_ASSEMBLIES   } from '../../modules/local/downloadassemblies'

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

    emit:
    stats           = FIND_ASSEMBLIES.out.stats        // channel: [ val(taxon), file(stats) ]
    versions        = ch_versions                      // channel: [ versions.yml ]
}
