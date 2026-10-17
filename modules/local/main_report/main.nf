process MAIN_REPORT {
    tag "$group_meta.id"
    label 'process_low'

    conda "conda-forge::quarto=1.6.41 bioconda::r-pathosurveilr=0.4.8"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e4/e487168aaa5f7b7a2dfabc1f308869fac1c466ab6bd34acd7a3715d927337c74/data':
        'community.wave.seqera.io/library/r-pathosurveilr_quarto:d4f39be8e8ae4734' }"

    input:
    tuple val(group_meta), file(inputs), path(template, stageAs: 'main_report_template')

    output:
    tuple val(group_meta), path("${prefix}.html"), emit: html
    tuple val(group_meta), path("${prefix}.pdf") , emit: pdf, optional: true
    path "versions.yml"               , emit: versions_main_report

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def tmpl = (group_meta.template ?: '').toString().trim()
    // A template may be given as an absolute path; use only its final segment for naming.
    // Interpolating the raw value made both --output-dir and the cp target multi-segment, and
    // quarto's nested output tree then failed the cp with "No such file or directory". Taking
    // the basename also keeps the bare-name case consistent: `pathsurveil_dashboard` and
    // /abs/path/pathsurveil_dashboard name the report identically. The published file is
    // ${prefix}.html, so the "report" keyword is carried by the directory name itself
    // (pathsurveil_report) and must not be stripped from the label.
    // Strip any trailing slashes first, so a harmless typo such as /abs/pathsurveil_report/ still
    // names the report after the directory. Note the basename cannot come from template.name:
    // the input is staged as 'main_report_template', so that would yield the staging name.
    def label = tmpl.replaceAll('/+$', '')
    label = label.contains('/') ? label.substring(label.lastIndexOf('/') + 1) : label
    if (label == '.' || label == '..') {
        throw new IllegalArgumentException("report_data template '${tmpl}' does not name a usable report directory")
    }
    prefix = task.ext.prefix ?: "${group_meta.id}${label ? '_' + label : ''}"
    """
    # Needed to avoid this issue: https://github.com/conda-forge/quarto-feedstock/issues/30
    if [[ -f /opt/conda/etc/conda/activate.d/quarto.sh ]]; then
        source /opt/conda/etc/conda/activate.d/quarto.sh
    fi
    # Tell quarto where to put cache so it does not try to put it where it does nmt have permissions
    export XDG_CACHE_HOME="\$(pwd)/cache"

    # Template is always a directory containing .qmd (dir-only per report_data spec)
    cp -r --dereference main_report_template main_report
    # Drop local development state that should never be staged: .quarto/ is a developer-machine
    # quarto project cache (xref/idx/freeze) and .gitignore is only meaningful in the source tree.
    rm -rf main_report/.quarto main_report/.gitignore

    # Render the report
    # NOTE: quarto resolves --output-dir relative to the project dir, so the rendered site
    # lands in main_report/${prefix}/ (and a website project nests it under _site/).
    quarto render main_report \\
        ${args} \\
        --output-dir "${prefix}" \\
        -P inputs:../${inputs}

    # Locate the rendered report page.
    # NOTE: deliberately no find/xargs/grep here. The task image ships no GNU findutils, so
    # find and xargs are absent -- xargs reported that as a bare 127 while a swallowed stderr
    # hid the real message. Globs are expanded by bash, so only coreutils are needed.
    for tool in ls head cp; do
        command -v \$tool >/dev/null 2>&1 || { echo "ERROR: required tool '\$tool' not found in the task image" >&2; exit 1; }
    done
    # NOTE: every use of ${prefix} below is quoted. A directory name may contain spaces
    # (e.g. "my templates"), and an unquoted expansion would word-split out_dir, silently
    # break the glob, and make cp target the wrong path.
    out_dir="main_report/${prefix}"
    if [[ ! -d \$out_dir ]]; then
        echo "ERROR: expected quarto output directory \$out_dir was not created" >&2
        exit 1
    fi
    # Largest page wins ('ls -S'): a website project's index.html can be a small redirect stub
    # sitting next to the substantive report page. nullglob plus plain globs keep this working
    # under the bash 3.2 shipped with macOS, so no globstar is used.
    # NOTE: each glob is double-quoted per component, because bash splits an unquoted variable
    # on IFS (spaces in a directory name would otherwise truncate the pattern) and expands the
    # glob against whatever split fragments resulted.
    shopt -s nullglob
    html_files=( "\$out_dir"/*.html "\$out_dir"/*/*.html "\$out_dir"/*/*/*.html )
    shopt -u nullglob
    if [[ \${#html_files[@]} -eq 0 ]]; then
        echo "ERROR: quarto produced no HTML in \$out_dir" >&2
        exit 1
    fi
    real_report=\$(ls -S "\${html_files[@]}" | head -n1)
    cp -L "\$real_report" "${prefix}.html"

    # Clean up
    rm -r main_report

    # Save version of quarto used
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto --version)
        r-PathoSurveilR: \$(Rscript -e "cat(as.character(packageVersion('dplyr')))")
    END_VERSIONS
    """
}
