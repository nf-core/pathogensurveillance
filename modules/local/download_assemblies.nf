process DOWNLOAD_ASSEMBLIES {                                                  
    tag "$id"                                                              
    label 'process_single'
    maxForks 1
    errorStrategy { return task.attempt > 3 ? 'ignore' : 'retry' }
    maxRetries 5                                                      
                                                                                
    conda "conda-forge::ncbi-datasets-cli=15.11.0"                                
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-datasets-cli:14.26.0 ':
        'docker.io/biocontainers/ncbi-datasets-cli:14.2.2_cv2' }"               
                                                                                
    input:                                                                      
    val id // There is no meta because we want to cache based only the ID                                                                
                                                                                
    output:                                                                     
    tuple val(id), path("${prefix}.zip"), emit: assembly
    tuple val(id), path("ncbi_dataset/data/${prefix}/${prefix}_*_genomic.fna"), emit: sequence
    tuple val(id), path("ncbi_dataset/data/${prefix}/${prefix}.gff"), emit: gff, optional: true
    tuple val(id), path("ncbi_dataset/data/${prefix}/${prefix}_cds.fna"), emit: cds, optional: true
    tuple val(id), path("ncbi_dataset/data/${prefix}/${prefix}.faa"), emit: protein, optional: true
                             
    path "versions.yml"                 , emit: versions                       
                                                                                
    when:                                                                       
    task.ext.when == null || task.ext.when                                      
                                                                                
    script:                                                                     
    prefix = task.ext.prefix ?: "${id}".replaceAll(' ', '_')                               
    def args = task.ext.args ?: ''                                              
    """
    # Download assemblies as zip archives
    datasets download genome accession $id --include gff3,rna,cds,protein,genome,seq-report --filename ${prefix}.zip
    
    # Unzip
    unzip ${prefix}.zip

    # Rename files with assembly name
    if [ -f ncbi_dataset/data/${prefix}/genomic.gff ]; then
        mv ncbi_dataset/data/${prefix}/genomic.gff ncbi_dataset/data/${prefix}/${prefix}.gff
    fi
    if [ -f ncbi_dataset/data/${prefix}/cds_from_genomic.fna ]; then
        mv ncbi_dataset/data/${prefix}/cds_from_genomic.fna ncbi_dataset/data/${prefix}/${prefix}_cds.fna
    fi
    if [ -f ncbi_dataset/data/${prefix}/protein.faa ]; then
        mv ncbi_dataset/data/${prefix}/protein.faa ncbi_dataset/data/${prefix}/${prefix}.faa
    fi
 
    cat <<-END_VERSIONS > versions.yml                                          
    "${task.process}":                                                          
        datasets: \$(datasets --version | sed -e "s/datasets version: //")                                           
    END_VERSIONS                                                                
    """                                                                         
}
