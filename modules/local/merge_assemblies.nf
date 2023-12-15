process MERGE_ASSEMBLIES {
    tag "All"
    label 'process_single'
                                                                                
    conda "conda-forge::coreutils=9.1"                                          
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :            
        'nf-core/ubuntu:20.04' }"                                                       

    input:
    path stats // multiple outputs from FIND_ASSEMBLIES

    output:
    path "merged_assembly_stats.tsv", emit: merged_stats

    when:
    task.ext.when == null || task.ext.when

    script:
    def first = stats[0]
    """
    # Save header from first path, adding a header for a new column for family
    head -n 1 ${first} | sed 's/\$/\\tFamily/' > merged_assembly_stats.tsv

    # Copy all data to this new file, except for headers, adding family column contents
    find ${stats} | xargs -I '{}' sh -c 'tail -n +2 {} | sed "s/\$/\\t\$(basename {} .tsv)/"' >> merged_assembly_stats.tsv
    """
}
