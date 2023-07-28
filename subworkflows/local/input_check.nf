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
        .set { reads_and_ref }
    
    emit:
    reads_and_ref                             // channel: [ val(meta), [ file(reads) ], val(ref_meta), file(ref), [ val(groups) ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_reads_ref_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.single_end = row.single_end.toBoolean()

    // check for presence of input files
    if (!file(row.reference).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Reference file does not exist!\n${row.reference}"
    }
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }

    // add path(s) of the fastq file(s) to the meta map
    def ref_meta = [:]
    ref_meta.id = row.reference_id ?: row.reference.replaceAll('/', '_')
    def output = []
    if (meta.single_end) {
        output = [ meta, [ file(row.fastq_1) ], ref_meta, file(row.reference), row.report_group.split(";") ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        output = [ meta, [ file(row.fastq_1), file(row.fastq_2) ], ref_meta, file(row.reference), row.report_group.split(";") ]
    }
    return output
}


