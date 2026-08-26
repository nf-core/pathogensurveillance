process FIND_ASSEMBLIES {
    tag "$taxon"
    label 'process_single'

    conda "conda-forge::ncbi-datasets-cli=18.31.0"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/ncbi-datasets-cli:18.31.0--73fbbae7aa16069b':
        'community.wave.seqera.io/library/ncbi-datasets-cli:18.31.0--f7fda2139b40106c' }"

    input:
    val taxon

    output:
    tuple val(taxon), path(output_path), emit: stats
    path "versions.yml"                , emit: versions_find_assemblies, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: taxon.replaceAll(' ', '_')
    args = task.ext.args ?: ''
    split_taxonomy = task.ext.split_taxonomy ?: false
    taxonomy_args = task.ext.taxonomy_args ?: ''
    output_path = "${prefix}.json"
    """
    # Ensure that CA certificates are in path when docker/singularity are used
    if [ -f /opt/conda/ssl/cacert.pem ]; then
        export SSL_CERT_FILE=/opt/conda/ssl/cacert.pem
    fi

    if [ "${split_taxonomy}" = "true" ]; then
        # enumerate one-level child taxa to query separately
        (datasets summary taxonomy taxon ${taxonomy_args} ${taxon.toLowerCase()} 1> children.json 2> >(tee error.txt >&2); echo \$? > exit_code.txt) || true
        children=\$(grep -oE '"children":[0-9, ]+' children.json | grep -oE '[0-9]+' || true)
        for child in ${taxon.toLowerCase()} \$children; do
            # NOTE: This command errors when a taxon is found but has no data rather than just outputing an empty file,
            #   so the below code forces it to not fail and then fails if any other error occur
            (datasets summary genome taxon ${args} "\$child" 1>> ${output_path} 2> >(tee error.txt >&2); echo \$? > exit_code.txt) || true
            if [[ \$(cat exit_code.txt) -ne 0 ]] && [ -s error.txt ] && ! grep -q 'no genome data is currently available for this taxon.' error.txt; then
                exit 1
            fi
        done
    else
        (datasets summary genome taxon ${args} ${taxon.toLowerCase()} 1> ${output_path} 2> >(tee error.txt >&2); echo \$? > exit_code.txt) || true
        if [[ \$(cat exit_code.txt) -ne 0 ]] && [ -s error.txt ] && ! grep -q 'no genome data is currently available for this taxon.' error.txt; then
            exit 1
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        datasets: \$(datasets --version | sed -e "s/datasets version: //")
    END_VERSIONS
    """
}
