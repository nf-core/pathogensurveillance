process FIND_ASSEMBLIES {
    tag "$taxon"
    label 'process_single'

    conda "bioconda::entrez-direct=16.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/entrez-direct:16.2--he881be0_1':
        'biocontainers/entrez-direct:16.2--he881be0_1' }"

    input:
    val taxon // There is no meta because we dont want to cache based only the taxon
    val allow_non_refseq
    val allow_partial_refs
    val allow_unannotated

    output:
    tuple val(taxon), path("${output_path}"), emit: stats
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${taxon}".replaceAll(' ', '_')
    def args = task.ext.args ?: ''
    output_path = "${prefix}"
    def filter_args = 'latest[PROP]'
    if (allow_partial_refs) {
        output_path = "${output_path}--partial"
    } else {
        filter_args = "${filter_args} AND full-genome-representation[PROP]"
    }
    if (allow_unannotated) {
        output_path = "${output_path}--unannotated"
    } else {
        filter_args = "${filter_args} AND has-annotation[PROP]"
    }
    if (allow_non_refseq) {
        output_path = "${output_path}--non_refseq"
    } else {
        filter_args = "${filter_args} NOT excluded-from-refseq[PROP]"
    }
    output_path = "${output_path}.tsv"
    """
    COLS="Id LastMajorReleaseAccession AssemblyName Taxid Organism SpeciesTaxid \\
          SpeciesName AssemblyType AssemblyStatus Coverage \\
          PartialGenomeRepresentation ReleaseLevel ReleaseType RefSeq_category \\
          SubmissionDate LastUpdateDate ContigN50 ScaffoldN50"

    echo \$COLS | sed 's/ /\\t/g' > ${output_path}
    esearch -db taxonomy -query "${taxon} OR ${taxon}[subtree]" | \\
        elink -target assembly | \\
        efilter -query "${filter_args}" | \\
        efetch -format docsum | \\
        xtract -pattern DocumentSummary -def 'NA' -element \$COLS >> \\
        ${output_path}

    if [[ \$(wc -l <${output_path}) -le 1 ]]; then
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esearch: \$(esearch -version)
    "${task.process}":
        elink: \$(elink -version)
    "${task.process}":
        efilter: \$(efilter -version)
    "${task.process}":
        efetch: \$(efetch -version)
    "${task.process}":
        xtract: \$(xtract -version)
    END_VERSIONS
    """
}
