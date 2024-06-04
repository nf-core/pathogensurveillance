//
// Check input samplesheets and convert to channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    sample_data_csv
    reference_data_csv

    main:

    SAMPLESHEET_CHECK ( sample_data_csv, reference_data_csv )

    sample_data = SAMPLESHEET_CHECK.out.sample_data
        .splitCsv ( header:true, sep:',', quote:'"' )
        .map { create_sample_metadata_channel(it) }
    reference_data = SAMPLESHEET_CHECK.out.reference_data
        .splitCsv ( header:true, sep:',', quote:'"' )
        .map { create_reference_metadata_channel(it) }

    sample_data = sample_data
        .transpose(by: 3)
        .combine(reference_data, by: 3)
        .map { ref_group_id, sample_meta, sample_paths, ncbi_accession, report_group_id, ref_meta, ref_path, ref_ncbi_accession ->
             [sample_meta, sample_paths, ncbi_accession, report_group_id, ref_meta]
        }
        .groupTuple(by: 0..3)

    emit:
    sample_data       // sammple_meta, [paths], ncbi_accession, report_meta, [ref_meta]
    reference_data    // ref_meta, path, ncbi_accession, ref_group_id
    sample_metadata_csv = SAMPLESHEET_CHECK.out.sample_data
    reference_metadata_csv = SAMPLESHEET_CHECK.out.reference_data
    versions = SAMPLESHEET_CHECK.out.versions
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


