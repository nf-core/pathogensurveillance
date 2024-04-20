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
        .map { create_reads_ref_channel(it) } // meta, [reads], reads_sra, ref_meta, reference, reference_refseq, groups
        .transpose ( by: 6 ) // Duplicate rows for each group when there are multiple groups per sample
        .map { it[0..5] + [[id: it[6]]] }
        .set { sample_data }

    emit:
    sample_data                               // meta, [reads], reads_sra, ref_meta, reference, reference_refseq, group
    csv      = SAMPLESHEET_CHECK.out.csv      // modified csv of metadata
    versions = SAMPLESHEET_CHECK.out.versions // versions.yml
}

def create_reads_ref_channel(LinkedHashMap row) {
    // check for presence of input files
    if (row.reference != '' && !file(row.reference).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Reference file does not exist.\n${row.reference}"
    }
    if (row.reads_1 != "" && !file(row.reads_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Reads 1 FastQ file does not exist.\n${row.reads_1}"
    }
    if (row.reads_2 != "" && !file(row.reads_2).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Reads 2 FastQ file does not exist.\n${row.reads_2}"
    }

    // Create sample metadata map
    def meta = [id: row.sample_id ?: null, single_end: row.reads_2 == '', reads_type: row.reads_type]

    // Create reference metadata map
    def ref_meta = [id: row.reference_id ?: null]

    // Format path(s) of the reads fastq file(s)
    def reads = null
    if (row.reads_1 != "") {
        reads = [file(row.reads_1)]
        if (row.reads_2 != "") {
            reads.add(file(row.reads_2))
        }
    }

    // Format paths to single files
    def sra = row.reads_sra ?: null
    def reference = row.reference ? file(row.reference): null
    def reference_refseq = row.reference_refseq ?: null
    def groups = row.report_group.split(";") as ArrayList

    return [meta, reads, sra, ref_meta, reference, reference_refseq, groups]
}


