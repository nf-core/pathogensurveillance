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
        .map { create_reads_ref_channel(it) } // meta, [shortread], nanopore, sra, ref_meta, reference, reference_refseq, groups
        .transpose ( by: 7 ) // Duplicate rows for each group when there are multiple groups per sample
        .map { it[0..6] + [[id: it[7]]] }
        .set { sample_data }

    emit:
    sample_data                               // meta, [shortread], nanopore, sra, ref_meta, reference, reference_refseq, group
    csv      = SAMPLESHEET_CHECK.out.csv      // modified csv of metadata
    versions = SAMPLESHEET_CHECK.out.versions // versions.yml
}

def create_reads_ref_channel(LinkedHashMap row) {
    // check for presence of input files
    if (row.reference != '' && !file(row.reference).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Reference file does not exist.\n${row.reference}"
    }
    if (row.shortread_1 != "" && !file(row.shortread_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Shortread 1 FastQ file does not exist.\n${row.shortread_1}"
    }
    if (row.shortread_2 != "" && !file(row.shortread_2).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Shortread 2 FastQ file does not exist.\n${row.shortread_2}"
    }
    if (row.nanopore != "" && !file(row.nanopore).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Nanopore FastQ file does not exist.\n${row.nanopore}"
    }

    // Create sample metadata map
    def meta = [id: row.sample_id ?: null, single_end: row.shortread_2 == '']

    // Create reference metadata map
    def ref_meta = [id: row.reference_id ?: null]
    
    // Format path(s) of the shortread fastq file(s) 
    def shortread = null
    if (row.shortread_1 != "") {
        shortread = [row.shortread_1]
        if (row.shortread_2 != "") {
            shortread.add(row.shortread_2)
        }
    }
    
    // Format paths to single files
    def nanopore = row.nanopore ?: null
    def sra = row.sra ?: null
    def reference = row.reference ?: null
    def reference_refseq = row.reference_refseq ?: null
    def groups = row.report_group.split(";") as ArrayList
    
    return [meta, shortread, nanopore, sra, ref_meta, reference, reference_refseq, groups]
}


