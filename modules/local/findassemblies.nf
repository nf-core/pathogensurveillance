process FIND_ASSEMBLIES {                                                  
    tag "$taxon"                                                              
    label 'process_single'
    maxForks 1                                                      
                                                                                
    conda "bioconda::entrez-direct=16.2"                                        
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/entrez-direct:16.2--he881be0_1':
        'quay.io/biocontainers/entrez-direct:16.2--he881be0_1' }"               
                                                                                
    input:                                                                      
    val taxon // There is no meta because we dont want to cache based only the taxon                                                                
                                                                                
    output:                                                                     
    tuple val(taxon), path("${prefix}.tsv"), emit: stats                                  
    path "versions.yml"                    , emit: versions                             
                                                                                
    when:                                                                       
    task.ext.when == null || task.ext.when                                      
                                                                                
    script:                                                                     
    prefix = task.ext.prefix ?: "${taxon}".replaceAll(' ', '_')                               
    def args = task.ext.args ?: ''                                              
    """
    COLS="Id LastMajorReleaseAccession AssemblyName Taxid Organism SpeciesTaxid \\
          SpeciesName AssemblyType AssemblyStatus Coverage \\
          PartialGenomeRepresentation ReleaseLevel ReleaseType RefSeq_category \\
          SubmissionDate LastUpdateDate ContigN50 ScaffoldN50"
    
    echo \$COLS | sed 's/ /\\t/g' > ${prefix}.tsv                                                                     
    esearch -db taxonomy -query "${taxon} OR ${taxon}[subtree]" | \\
        elink -target assembly | \\
        efilter -query "latest [PROP] NOT partial-genome-representation [PROP]" | \\
        efetch -format docsum | \\
        xtract -pattern DocumentSummary -def 'NA' -element \$COLS >> \\
        ${prefix}.tsv
    
    cat <<-END_VERSIONS > versions.yml                                          
    "${task.process}":                                                          
        esearch: \$(esearch -version)                                           
    "${task.process}":                                                          
        elink: \$(elink -version)                                           
    "${task.process}":                                                          
        efilter: \$(efilter -version)                                           
    "${task.process}":                                                          
        efetch: \$(efetch -version)                                           
    "${task.process}":                                                          
        xtract: \$(xtract -version)                                           
    END_VERSIONS                                                                
    """                                                                         
}
