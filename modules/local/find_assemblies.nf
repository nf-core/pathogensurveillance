process FIND_ASSEMBLIES {
    tag "$taxon"
    label 'process_single'

    conda "conda-forge::ncbi-datasets-cli=16.35.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://logan-blair/collection/ncbi-datasets-cli:latest':
        'docker.io/staphb/ncbi-datasets:16.35.0' }"

    input:
    val taxon // There is no meta because we dont want to cache based only the taxon

    output:
    tuple val(taxon), path(output_path), emit: stats
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: taxon.replaceAll(' ', '_')
    args = task.ext.args ?: ''
    output_path = "${prefix}.json"
    """
    datasets summary genome taxon ${args} ${taxon.toLowerCase()} > ${output_path}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        datasets: \$(datasets --version | sed -e "s/datasets version: //")
    END_VERSIONS
    """
}
