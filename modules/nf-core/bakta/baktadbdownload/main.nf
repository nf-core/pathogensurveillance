process BAKTA_BAKTADBDOWNLOAD {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bakta:1.11.4--pyhdfd78af_0'
        : 'biocontainers/bakta:1.11.4--pyhdfd78af_0'}"

    output:
    path "db*"              , emit: db
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # TODO: Remove this workaround when bakta adds a proper User-Agent header on next release
    python3 -c "
    # Patch requests first
    import requests
    _original_request = requests.Session.request
    def patched_request(self, method, url, **kwargs):
        headers = kwargs.get('headers', {})
        if headers is None:
            headers = {}
        if 'User-Agent' not in headers:
            headers['User-Agent'] = 'bakta/1.11.4'
        kwargs['headers'] = headers
        return _original_request(self, method, url, **kwargs)
    requests.Session.request = patched_request
    
    # Import and run bakta
    import sys
    from bakta.db import main as bakta_main
    
    # Set up command line arguments
    sys.argv = ['bakta_db', 'download'] + '''$args'''.split()
    
    # Run it
    bakta_main()
    "

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$(echo \$(bakta_db --version) 2>&1 | cut -f '2' -d ' ')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    """
    echo "bakta_db \\
        download \\
        $args"

    mkdir db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$(echo \$(bakta_db --version) 2>&1 | cut -f '2' -d ' ')
    END_VERSIONS
    """
}
