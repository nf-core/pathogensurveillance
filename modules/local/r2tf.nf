process R2TF {
    tag "$ref_meta.id"
    label 'process_single'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(ref_meta), path(single_copy_fna), path(single_copy_faa)

    output:
    tuple val(ref_meta), path("${outdir}"), emit: output

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${ref_meta.id}"
    outdir = "${prefix}_r2tf_out"
    """
    # Defining codex
    CODEX=\$(echo "${ref_meta.name}" | awk '{print toupper(substr(\$1,1,3))toupper(substr(\$2,1,2))}')
    mkdir ${outdir}

    fna=\$(grep '^>' $single_copy_fna)

    while IFS= read -r fastafile; do
        echo "Processing file \$fastafile"

        # Extracting file name
        fasta=\$(echo \$fastafile | awk -F ':' '{print \$1}')
        buscoid=\$(basename \$fasta .fna)
        geneid="\${CODEX}-\${buscoid}"
    	header=\$(grep '>' \$fasta)

        # Genearting a unique IDentifier
        #NUM=\$(echo \$fasta | cut -d 'a' -f1)
        NUM=\$(echo \$header | awk -F '>' '{print \$2}' | awk -F ':' '{print \$1}')
	    echo "NUM \$NUM"
        ID="\$CODEX \$NUM"
        echo "Identifier: \$ID"

        # Replace the headers in the FASTA file
        echo "Processing file fasta \$fasta with header \$header"
        sed "s/\$header/>\$geneid \\| \$buscoid \\| \$geneid \\| [${ref_meta.name} ${ref_meta.id}]/" "\$fasta" | sed -e 's/\\(.*\\)/\\U\\1/'  >> ${outdir}/${prefix}_ntformatted.fa
        #printf "TAG\\n\\n" >> ${outdir}/${prefix}_ntformatted.fa
    done <<< \$fna


    faa=\$(grep '^>' $single_copy_faa)

    while IFS= read -r fastafile; do
        #echo "Processing file \$fastafile"

        # Extracting file name
        fasta=\$(echo \$fastafile | awk -F ':' '{print \$1}')
	    buscoid=\$(basename \$fasta .faa)
        geneid="\${CODEX}-\${buscoid}"
        header=\$(grep '>' \$fasta)

        # Genearting a unique IDentifier
        #NUM=\$(echo \$fasta | cut -d 'a' -f1)
        NUM=\$(echo \$header | awk -F '>' '{print \$2}' | awk -F ':' '{print \$1}')
    	echo "NUM \$NUM"
        ID="\$CODEX \$NUM"
        echo "Identifier: \$ID"

        # Replace the headers in the FASTA file
        echo "Processing file fasta \$fasta with header \$header"
        sed "s/\$header/>\$geneid \\| \$buscoid \\| \$geneid \\| [${ref_meta.name} ${ref_meta.id}]/" "\$fasta" | sed -e 's/\\(.*\\)/\\U\\1/' >> ${outdir}/${prefix}_aaformatted.fa
    done <<< \$faa
    """
}
