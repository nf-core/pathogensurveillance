include { FIND_ASSEMBLIES       } from '../../modules/local/findassemblies'
include { PICK_ASSEMBLIES       } from '../../modules/local/pickassemblies'

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

    PICK_ASSEMBLIES (
        ch_families                                                       
        .join(ch_genera)                                                        
        .join(ch_species)                                                       
        .splitText( elem: 1 )                                                   
        .transpose( by: 1 )                                                     
        .join(FIND_ASSEMBLIES.out.stats, by: 1)
    )


    emit:
    stats           = FIND_ASSEMBLIES.out.stats        // channel: [ val(taxon), file(stats) ]
    versions        = ch_versions                      // channel: [ versions.yml ]
}
