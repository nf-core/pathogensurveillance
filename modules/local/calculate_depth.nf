process CALCULATE_DEPTH {                                                          
    tag "$meta.id"                                                              
    label 'process_single'                                                      
                                                                                
    conda "conda-forge::coreutils=9.1"                                          
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :            
        'nf-core/ubuntu:20.04' }"                                                       
                                                                                
    input:                                                                      
    tuple val(meta), path(fastqs), path(ref)                                     
                                                                                
    output:                                                                     
    tuple val(meta), env(DEPTH), emit: depth                     
                                                                                
    when:                                                                       
    task.ext.when == null || task.ext.when                                      
                                                                                
    script:                                                                     
    prefix = task.ext.prefix ?: "${meta.id}"                                    
    def args = task.ext.args ?: ''                                              
    """                                                                         
    REF_WC=\$(zgrep -v '^>' ${ref} | wc -m)                                     
    READ_COUNT=\$(zgrep -c '@' ${fastqs[0]})                                    
    READ_LEN=\$(zgrep -m 1 '^[^@+#]' ${fastqs[0]} | wc -m)                      
    DEPTH=\$((\$READ_COUNT * \$READ_LEN / \$REF_WC))                            
    """                                                                         
}                                                                               

