//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
    SAMPLESHEET_CHECK.out
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_reads_ref_channel(it) }
        .transpose ( by: 4 ) // Duplicate rows for each group when there are multiple groups per sample
        .map { it[0..3] + [[id: it[4]]] }
        .set { sample_data }

    emit:
    sample_data                             // channel: [ val(meta), [ file(reads) ], val(ref_meta), file(ref), val(groups) ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_reads_ref_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.single_end = row.single_end.toBoolean()

    // check for presence of input files
    if (row.reference != '' && !file(row.reference).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Reference file does not exist!\n${row.reference}"
    }
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }

    // add path(s) of the fastq file(s) to the meta map
    def ref_meta = [id: null]
    def output = []
    def groups = row.report_group.split(";") as ArrayList
    def reference = null
    if (row.reference != '') {
        reference = file(row.reference)
        ref_meta.id = row.reference_id ?: row.reference.replaceAll('/', '_')
    }
    if (meta.single_end) {
        output = [ meta, [ file(row.fastq_1) ], ref_meta, reference, groups ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        output = [ meta, [ file(row.fastq_1), file(row.fastq_2) ], ref_meta, reference, groups ]
    }
    return output
}


