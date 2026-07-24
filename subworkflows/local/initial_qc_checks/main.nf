include { FASTQC        } from '../../../modules/nf-core/fastqc/main'
include { NANOPLOT      } from '../../../modules/nf-core/nanoplot/main'

workflow INITIAL_QC_CHECKS {

    take:
    sample_data

    main:
    versions = channel.empty()
    messages = channel.empty()

    // Run FastQC
    shortreads = sample_data
        .filter { sample_meta -> sample_meta.sequence_type == "illumina" || sample_meta.sequence_type == "bgiseq" }
        .map { sample_meta -> [[id: sample_meta.sample_id], sample_meta.paths] }
        .unique()
    FASTQC ( shortreads )

    // Run Nanoplot
    nanopore_reads = sample_data
        .filter { sample_meta -> sample_meta.sequence_type == "nanopore" || sample_meta.sequence_type == "pacbio" }
        .map { sample_meta -> [[id: sample_meta.sample_id], sample_meta.paths] }
        .unique()
    NANOPLOT ( nanopore_reads )

    emit:
    fastqc_zip    = FASTQC.out.zip
    nanoplot_txt  = NANOPLOT.out.txt
    versions      = versions         // versions
    messages      = messages         // meta, group_meta, ref_meta, workflow, level, message
}
