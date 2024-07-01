process DOWNLOAD_ASSEMBLIES {
    tag "${ref_meta.id}"
    label 'process_single'
    maxForks 3
    errorStrategy { return task.attempt > 2 ? 'ignore' : 'retry' }
    maxRetries 6
    maxErrors 10

    conda "conda-forge::ncbi-datasets-cli=15.11.0 conda-forge::unzip=6.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-datasets-cli:14.26.0 ':
        'zachary-foster/ncbi-datasets-cli:16.0.0' }"

    input:
    tuple val(ref_meta), val(id)

    output:
    tuple val(ref_meta), path("${prefix}.zip"), emit: assembly
    tuple val(ref_meta), path("ncbi_dataset/data/${id}/${prefix}_genomic.fna"), emit: sequence
    tuple val(ref_meta), path("ncbi_dataset/data/${id}/${prefix}.gff"), emit: gff, optional: true
    tuple val(ref_meta), path("ncbi_dataset/data/${id}/${prefix}_cds.fna"), emit: cds, optional: true
    tuple val(ref_meta), path("ncbi_dataset/data/${id}/${prefix}.faa"), emit: protein, optional: true

    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${ref_meta.id}"
    def args = task.ext.args ?: ''
    """
    # Test that unzip is in the docker image
    # unzip --version

    # Download assemblies as zip archives
    datasets download genome accession $id --include gff3,rna,cds,protein,genome,seq-report --filename ${prefix}.zip

    # Unzip
    unzip ${prefix}.zip

    # Rename files with assembly name
    if [ -f ncbi_dataset/data/${id}/${id}_*_genomic.fna ]; then
        mv ncbi_dataset/data/${id}/${id}_*_genomic.fna ncbi_dataset/data/${id}/${prefix}_genomic.fna
    fi
    if [ -f ncbi_dataset/data/${id}/genomic.gff ]; then
        mv ncbi_dataset/data/${id}/genomic.gff ncbi_dataset/data/${id}/${prefix}.gff
    fi
    if [ -f ncbi_dataset/data/${id}/cds_from_genomic.fna ]; then
        mv ncbi_dataset/data/${id}/cds_from_genomic.fna ncbi_dataset/data/${id}/${prefix}_cds.fna
    fi
    if [ -f ncbi_dataset/data/${id}/protein.faa ]; then
        mv ncbi_dataset/data/${id}/protein.faa ncbi_dataset/data/${id}/${prefix}.faa
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        datasets: \$(datasets --version | sed -e "s/datasets version: //")
    END_VERSIONS
    """
}
