process REFORMAT_PIRATE_RESULTS {
    tag "$ref_meta.id"
    label 'process_single'
                                                                                
    conda "bioconda::pirate=1.0.5 bioconda::perl-bioperl=1.7.8"                 
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pirate:1.0.4--hdfd78af_2' :
        'quay.io/biocontainers/pirate:1.0.5--hdfd78af_0' }"                     

    input:
    tuple val(ref_meta), path(pirate_results)

    output:
    tuple val(ref_meta), path("${prefix}.tsv")                , emit: gene_fam
    tuple val(ref_meta), path("${prefix}_genePA.tsv")         , emit: gene_fam_pa
    tuple val(ref_meta), path("${prefix}_genePAparalogs.tsv") , emit: gene_fam_para
    tuple val(ref_meta), path("${prefix}_roary.tsv")          , emit: gene_fam_roary
    path "versions.yml"                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${ref_meta.id}"
    """
    # Reformat comments so that subsample_outputs.pl understands them
    mkdir reformatted_gffs
    cp modified_gffs/*.gff reformatted_gffs/
    sed -i 's/^# /## /' reformatted_gffs/*.gff

    # rename with original locus tag from input files
    subsample_outputs.pl -i PIRATE.gene_families.ordered.tsv -g reformatted_gffs/ -o ${prefix}.tsv --field "prev_ID" --feature "CDS"

    # gene/allele presence-absence
    PIRATE_to_Rtab.pl -i ${prefix}.tsv -o ${prefix}_genePA.tsv 

    # paralog presence-absence (duplications = d, fission/fusions = ff)
    paralogs_to_Rtab.pl -i ${prefix}.tsv -o ${prefix}_genePAparalogs.tsv --type d 

    #create roary format output
    PIRATE_to_roary.pl -i ${prefix}.tsv -o ${prefix}_roary.tsv

    cat <<-END_VERSIONS > versions.yml                                          
    "${task.process}":                                                          
        pirate: \$( echo \$( PIRATE --version 2>&1) | sed 's/PIRATE //' )       
    END_VERSIONS                                                                 
    """
}
