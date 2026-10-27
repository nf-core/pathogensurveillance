process MAKE_GFF_WITH_FASTA {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data' :
        'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), path(sequence), path(gff)

    output:
    tuple val(meta), path("${prefix}.gff"), emit: gff
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Copy gff info, removing the last "###" line
    if [[ $gff == *.gz ]]; then
        gunzip -c ${gff} | head -n -1 > ${prefix}_with_ref.gff
    else
        head -n -1 ${gff} > ${prefix}_with_ref.gff
    fi

    # Add FASTA section header
    echo "##FASTA" >> ${prefix}_with_ref.gff

    # Add FASTA info, replacing headers with just ID
    if [[ $sequence == *.gz ]]; then
        gunzip -c ${sequence} | sed -E 's/^>([a-zA-Z0-9_.]+) +.*\$/>\\1/g' >> ${prefix}_with_ref.gff
    else
        sed -E 's/^>([a-zA-Z0-9_.]+) +.*\$/>\\1/g' ${sequence} >> ${prefix}_with_ref.gff
    fi

    # Rename output file to be just the sample ID and make sure input file does not have same name
    mv ${gff} input_${gff}
    mv ${prefix}_with_ref.gff ${prefix}.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version | head -n 1 | sed 's/sed (GNU sed) //')
    END_VERSIONS
    """
}
