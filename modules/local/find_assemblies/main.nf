process FIND_ASSEMBLIES {
    tag "$taxon"
    label 'process_single'

    conda "conda-forge::ncbi-datasets-cli=18.7.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fc/fc3b4f67d6b34c43dba4fd417994f88e2b5042c0d68a9e53d59de7174c90f119/data':
        'community.wave.seqera.io/library/ncbi-datasets-cli:18.7.0--03372988f820bf79' }"

    input:
    val taxon

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
    # NOTE: This command errors when a taxon is found but has no data rather than just outputing an empty file,
    #   so the below code forces it to not fail and then fails if any other error occur
    (datasets summary genome taxon ${args} ${taxon.toLowerCase()} 1> ${output_path} 2> >(tee error.txt >&2); echo \$? > exit_code.txt) || true
    if [[ \$(cat exit_code.txt) -ne 0 ]] && [ -s error.txt ] && ! grep -q 'no genome data is currently available for this taxon.' error.txt; then
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        datasets: \$(datasets --version | sed -e "s/datasets version: //")
    END_VERSIONS
    """
}
