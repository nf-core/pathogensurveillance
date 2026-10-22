process RENDER_REPORT {
    tag "$group_meta.id"
    label 'process_low'

    conda "conda-forge::quarto=1.6.41 bioconda::r-pathosurveilr=0.4.8 conda-forge::latexmk=4.88 conda-forge::tectonic=0.17.0"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/15/159334bee8e9c74c1c964125f9d586f97481457a63eab737c86894224c7c4c1b/data':
        'community.wave.seqera.io/library/r-pathosurveilr_latexmk_quarto_tectonic:d274308f89b69fd3' }"

    input:
    tuple val(group_meta), file(inputs), path(template, stageAs: 'render_report_template')

    output:
    tuple val(group_meta), path("${prefix}*"), emit: report, optional: true
    path "versions.yml", emit: versions_render_report

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def tmpl = (group_meta.template ?: '').toString().trim()
    def render_target = (group_meta.render_target ?: '').toString().trim()
    def label = tmpl.replaceAll('/+$', '').tokenize('/').last()
    prefix = task.ext.prefix ?: "${group_meta.id}${label ? '_' + label : ''}"
    def target_extension = render_target ? render_target.substring(render_target.lastIndexOf('.')) : ''
    def target_shell = "'" + render_target.replace("'", "'\"'\"'") + "'"
    """

    # Needed to avoid this issue: https://github.com/conda-forge/quarto-feedstock/issues/30
    if [[ -f /opt/conda/etc/conda/activate.d/quarto.sh ]]; then
        source /opt/conda/etc/conda/activate.d/quarto.sh
    fi

    # Tell quarto where to put cache so it does not try to put it where it does nmt have permissions
    export XDG_CACHE_HOME="\$(pwd)/cache"

    # Render a private copy; the staged template may be a symlink to user data.
    cp -r --dereference render_report_template report_template
    rm -rf report_template/.quarto

    quarto render report_template \\
        ${args} \\
        --output-dir "${prefix}" \\
        -P inputs:../${inputs}

    # Locate rendered output
    site_dir="report_template/${prefix}"
    if [[ -d "\$site_dir/_site" ]]; then
        site_dir="\$site_dir/_site"
    fi

    # Publish single target or full directory
    if [[ -n "${render_target}" ]]; then
        target_file="\$site_dir"/${target_shell}
        if [[ -f "\$target_file" ]]; then
            cp -- "\$target_file" "${prefix}${target_extension}"
        fi
    else
        if [[ -d "\$site_dir" ]]; then
            mkdir -p "${prefix}"
            cp -R "\$site_dir/." "${prefix}/"
        fi
    fi

    # Clean up
    rm -r report_template

    # Save version of quarto used
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto --version)
        r-PathoSurveilR: \$(Rscript -e "cat(as.character(packageVersion('PathoSurveilR')))")
    END_VERSIONS
    """
}
