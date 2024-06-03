//
// Check input samplesheets and convert to channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    
    SAMPLESHEET_CHECK ( samplesheet )
    
    sample_metadata = SAMPLESHEET_CHECK.out.sample_metadata
        .splitCsv ( header:true, sep:',', quote:'"' )
        .map { create_sample_metadata_channel(it) }
    reference_metadata = SAMPLESHEET_CHECK.out.reference_metadata
        .splitCsv ( header:true, sep:',', quote:'"' )
        .map { create_reference_metadata_channel(it) }

    emit:
    sample_metadata                           // sammple_meta, [paths], ncbi_accession, [ref_group_ids], report_group_id
    reference_metadata                        // ref_meta, path, ncbi_accession, ref_group_id
    csv      = SAMPLESHEET_CHECK.out.csv      // modified csv of metadata
    versions = SAMPLESHEET_CHECK.out.versions // versions.yml
}

def create_sample_metadata_channel(LinkedHashMap row) {
    if (row.path != "" && !file(row.path).exists()) {
        exit 1, "ERROR: Please check the sample metadata CSV. The file specified by 'path` does not exist.\n${row.path}"
    }
    if (row.path_2 != "" && !file(row.path_2).exists()) {
        exit 1, "ERROR: Please check the sample metadata CSV. The file specified by 'path_2' does not exist.\n${row.path_2}"
    }
    def sammple_meta = [
        id: row.sample_id,
        name: row.name,
        description: row.description,
        sequence_type: row.sequence_type,
        ploidy: row.ploidy,
        single_end: row.path_2 == ''
    ]
    def paths = null
    if (row.path != "") {
        paths = [file(row.path)]
        if (row.path_2 != "") {
            paths.add(file(row.path_2))
        }
    }
    def ref_group_ids = row.ref_group_ids.split(";") as ArrayList
    def ncbi_accession = row.ncbi_accession ?: null
    def report_group_id = [id: row.report_group_ids]
    return [sammple_meta, paths, ncbi_accession, ref_group_ids, report_group_id]
}

def create_reference_metadata_channel(LinkedHashMap row) {
    if (row.ref_path != "" && !file(row.ref_path).exists()) {
        exit 1, "ERROR: Please check the reference metadata CSV. The file specified by 'ref_path` does not exist.\n${row.ref_path}"
    }
    def ref_meta = [
        id: row.ref_id,
        name: row.ref_name,
        description: row.ref_description,
        primary_usage: row.ref_primary_usage,
        contextual_usage: row.ref_contextual_usage
    ]
    def path = row.ref_path
    def ref_group_id = row.ref_group_ids
    def ncbi_accession = row.ref_ncbi_accession ?: null
    return [ref_meta, path, ncbi_accession, ref_group_id]
}


