//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
    SAMPLESHEET_CHECK.out.csv
        .splitCsv ( header:true, sep:',', quote:'"')
        .map { create_reads_ref_channel(it) }
        .transpose ( by: 6 ) // Duplicate rows for each group when there are multiple groups per sample
        .map { it[0..5] + [[id: it[6]]] }
        .set { sample_data }

    emit:
    sample_data                               // val(meta), [ file(reads) ], val(sra), val(ref_meta), file(ref), val(ref_acc), val(groups)
    csv      = SAMPLESHEET_CHECK.out.csv      // modified csv of metadata
    versions = SAMPLESHEET_CHECK.out.versions // versions.yml
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ], sra ]
def create_reads_ref_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.single_end = row.single_end.toBoolean()

    // check for presence of input files
    if (row.reference != '' && !file(row.reference).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Reference file does not exist!\n${row.reference}"
    }
    if (row.fastq_1 != "" && !file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (row.reference != "" && row.reference_accession != "") {
        exit 1, "ERROR: Both a reference path (URL or local file) and a reference accession ID (e.g. from RefSeq) cannot be defined."
    }

    // add path(s) of the fastq file(s) to the meta map
    def ref_meta = [id: null]
    def output = []
    def groups = row.report_group.split(";") as ArrayList
    def reference = null
    def sra = row.sra
    def reference_accession = row.reference_accession ?: ''
    if (row.reference != '') {
        reference = file(row.reference)
        ref_meta.id = row.reference_id ?: row.reference.replaceAll('[\\/:*?"<>| .]', '_')
    }
    if (reference_accession != '') {
        ref_meta.id = row.reference_id ?: reference_accession.replaceAll('[\\/:*?"<>| .]', '_')
    }
    if (sra != '') {
        output = [ meta, [], sra, ref_meta, reference, reference_accession, groups ]
    } else if (meta.single_end) {
        output = [ meta, [ file(row.fastq_1) ], null, ref_meta, reference, reference_accession, groups ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        output = [ meta, [ file(row.fastq_1), file(row.fastq_2) ], null, ref_meta, reference, reference_accession, groups ]
    }
    return output
}


