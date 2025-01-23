include { FASTQC        } from '../../modules/nf-core/fastqc/main'
include { NANOPLOT      } from '../../modules/nf-core/nanoplot/main'

workflow INITIAL_QC_CHECKS {

    take:
    sample_data

    main:
    versions = Channel.empty()
    messages = Channel.empty()

    // Run FastQC
    shortreads = sample_data
        .filter { it.sequence_type == "illumina" || it.sequence_type == "bgiseq" }
        .map { [[id: it.sample_id], it.paths] }
        .unique()
    FASTQC ( shortreads )
    versions = versions.mix(FASTQC.out.versions)

    // Run Nanoplot
    nanopore_reads = sample_data
        .filter { it.sequence_type == "nanopore" || it.sequence_type == "pacbio" }
        .map { [[id: it.sample_id], it.paths] }
        .unique()
    NANOPLOT ( nanopore_reads )
    versions = versions.mix(NANOPLOT.out.versions)

    emit:
    fastqc_zip    = FASTQC.out.zip
    nanoplot_txt  = NANOPLOT.out.txt
    versions      = versions         // versions
    messages      = messages         // meta, group_meta, ref_meta, workflow, level, message
}
