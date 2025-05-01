process DOWNLOAD_ASSEMBLIES {
    tag "${ref_meta.id}"
    label 'process_single'

    conda "conda-forge::ncbi-datasets-cli=15.11.0 bioconda::samtools=1.18 conda-forge::unzip=6.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/nf-core/ncbi-datasets-cli:0.1':
        'quay.io/nf-core/ncbi-datasets-cli:0.1' }"

    input:
    tuple val(ref_meta), val(id)

    output:
    //tuple val(ref_meta), path("${prefix}.zip"), emit: assembly
    tuple val(ref_meta), path("${prefix}.fasta.gz"), emit: sequence
    tuple val(ref_meta), path("${prefix}.gff.gz"), emit: gff, optional: true
    //tuple val(ref_meta), path("${prefix}_cds.fna"), emit: cds, optional: true
    //tuple val(ref_meta), path("${prefix}.faa"), emit: protein, optional: true

    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${ref_meta.id}"
    def args = task.ext.args ?: ''
    """
    # Download assemblies as zip archives
    datasets download genome accession $id --include gff3,genome --filename ${prefix}.zip

    # Unzip
    unzip ${prefix}.zip

    # Rename files with assembly name
    if [ -f ncbi_dataset/data/${id}/${id}_*_genomic.fna ]; then
        mv ncbi_dataset/data/${id}/${id}_*_genomic.fna ${prefix}.fasta
    fi
    if [ -f ncbi_dataset/data/${id}/genomic.gff ]; then
        mv ncbi_dataset/data/${id}/genomic.gff ${prefix}.gff
    fi
    #if [ -f ncbi_dataset/data/${id}/cds_from_genomic.fna ]; then
    #    mv ncbi_dataset/data/${id}/cds_from_genomic.fna ${prefix}_cds.fna
    #fi
    #if [ -f ncbi_dataset/data/${id}/protein.faa ]; then
    #    mv ncbi_dataset/data/${id}/protein.faa ${prefix}.faa
    #fi

    # Delete unneeded files
    rm ${prefix}.zip
    rm -rf ncbi_dataset
    rm README.md md5sum.txt

    # Compress output files
    bgzip ${prefix}.fasta
    if [ -f ${prefix}.gff ]; then
        bgzip ${prefix}.gff
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        datasets: \$(datasets --version | sed -e "s/datasets version: //")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        unzip: \$(unzip -v | head -n 1 | awk '{print \$2}')
    END_VERSIONS
    """
}
