process R2TF {
    tag "$ref_meta.id"
    label 'process_single'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(ref_meta), path(busco_dir)

    output:
    tuple val(ref_meta), path("${outdir}"), emit: output

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${ref_meta.id}"
    outdir = "${prefix}_r2tf_out"
    """
    # Defining codex
    CODEX=${ref_meta.id}

    mkdir ${outdir}

    fna=\$(grep '^>' ${busco_dir}/\${CODEX}_genomic.fna/run_eukaryota_odb10/busco_sequences/single_copy_busco_sequences/*.fna)

    while IFS= read -r fastafile; do
        echo "Processing file \$fastafile"

        # Extracting file name
        fasta=\$(echo \$fastafile | awk -F ':' '{print \$1}')
        buscoid=\$(basename \$fasta .fna)
	header=\$(grep '>' \$fasta)
            
        # Genearting a unique IDentifier
        #NUM=\$(echo \$fasta | cut -d 'a' -f1)
        NUM=\$(echo \$header | awk -F '>' '{print \$2}' | awk -F ':' '{print \$1}')
	echo "NUM \$NUM" 
        ID="\$CODEX \$NUM"
        echo "Identifier: \$ID"

        # Replace the headers in the FASTA file
        echo "Processing file fasta \$fasta with header \$header"
        sed "s/\$header/>\$buscoid \\| [\$CODEX]/" "\$fasta" >> ${outdir}/\${CODEX}_ntformatted.fa
    done <<< \$fna


    faa=\$(grep '^>' ${busco_dir}/\${CODEX}_genomic.fna/run_eukaryota_odb10/busco_sequences/single_copy_busco_sequences/*.faa) 

    while IFS= read -r fastafile; do
        #echo "Processing file \$fastafile"

        # Extracting file name
        fasta=\$(echo \$fastafile | awk -F ':' '{print \$1}')
	buscoid=\$(basename \$fasta .faa)
        header=\$(grep '>' \$fasta)

        # Genearting a unique IDentifier
        #NUM=\$(echo \$fasta | cut -d 'a' -f1)
        NUM=\$(echo \$header | awk -F '>' '{print \$2}' | awk -F ':' '{print \$1}')
	echo "NUM \$NUM"
        ID="\$CODEX \$NUM"
        echo "Identifier: \$ID"

        # Replace the headers in the FASTA file
        echo "Processing file fasta \$fasta with header \$header"
        sed "s/\$header/>\$buscoid \\| [\$CODEX]/" "\$fasta" >> ${outdir}/\${CODEX}_aaformatted.fa
    done <<< \$faa
    """
}
