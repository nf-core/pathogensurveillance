//
// Check input samplesheets and convert to channels
//
include { BBMAP_SENDSKETCH  } from '../../modules/local/sendsketch'
include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'
include { SRATOOLS_FASTERQDUMP        } from '../../modules/local/fasterqdump'
include { INITIAL_CLASSIFICATION } from '../../modules/local/initial_classification'
include { DOWNLOAD_ASSEMBLIES                       } from '../../modules/local/download_assemblies'
include { FIND_ASSEMBLIES                           } from '../../modules/local/find_assemblies'
include { PICK_ASSEMBLIES                           } from '../../modules/local/pick_assemblies'

workflow PREPARE_INPUT {
    take:
    sample_data_csv
    reference_data_csv

    main:

    // Initalize channel to accumulate information about software versions used
    versions = Channel.empty()

    // Parse input CSVs
    SAMPLESHEET_CHECK ( sample_data_csv, reference_data_csv )
    sample_data = SAMPLESHEET_CHECK.out.sample_data
        .splitCsv ( header:true, sep:',', quote:'"' )
        .map { create_sample_metadata_channel(it) }
    reference_data = SAMPLESHEET_CHECK.out.reference_data
        .splitCsv ( header:true, sep:',', quote:'"' )
        .map { create_reference_metadata_channel(it) }

    // Download FASTQ files if an NCBI accession is provided
    ncbi_acc = sample_data
        .filter { sample_meta, sample_paths, ncbi_accession, ref_group_ids, report_meta ->
            ncbi_accession != null
        }
        .map { sample_meta, sample_paths, ncbi_accession, ref_group_ids, report_meta ->
            [sample_meta, ncbi_accession]
        }
        .unique()
    SRATOOLS_FASTERQDUMP ( ncbi_acc )
    sample_data = sample_data
        .combine(SRATOOLS_FASTERQDUMP.out.reads, by: 0)
        .map { sample_meta, sample_paths, ncbi_accession, ref_group_ids, report_meta, downloaded ->
            [sample_meta, downloaded ?: sample_paths, ncbi_accession, ref_group_ids, report_meta]
        }

    // Look up approximate taxonomic classifications
    reads = sample_data
        .map { sample_meta, sample_paths, ncbi_accession, ref_group_ids, report_meta ->
            [sample_meta, sample_paths]
        }
        .unique()
    BBMAP_SENDSKETCH ( reads )
    versions = versions.mix(BBMAP_SENDSKETCH.out.versions.toSortedList().map{it[0]})

    // Parse results of sendsketch to get list of taxa to download references for
    INITIAL_CLASSIFICATION ( BBMAP_SENDSKETCH.out.hits )
    versions = versions.mix(INITIAL_CLASSIFICATION.out.versions.toSortedList().map{it[0]})

    // Get list of families for all samples without exclusive references defined by the user
    all_families = INITIAL_CLASSIFICATION.out.families
        .map { sample_meta, families ->
            [families]
        }
        .splitText()
        .flatten()
        .map { families -> families.replace('\n', '') }
        .unique()

    // Download RefSeq metadata for all assemblies for every family found by the initial identification
    FIND_ASSEMBLIES ( all_families )
    versions = versions.mix(FIND_ASSEMBLIES.out.versions.toSortedList().map{it[0]})

    // Choose reference sequences to provide context for each sample
    PICK_ASSEMBLIES (
        INITIAL_CLASSIFICATION.out.families
            .join(INITIAL_CLASSIFICATION.out.genera)
            .join(INITIAL_CLASSIFICATION.out.species),
        FIND_ASSEMBLIES.out.stats
            .map { taxon, assembly_meta -> assembly_meta }
            .toSortedList(),
        params.n_ref_strains,
        params.n_ref_species,
        params.n_ref_genera
    )

    //new_reference_data = PICK_ASSEMBLIES.out.id_list
    //    //.map {sample_meta, new_ref_csv -> new_ref_csv}
    //    .splitText(elem: 1)
    //    .map { sample_meta, new_refs -> [sample_meta, new_ref_csv.replace('\n', '')] }
    //    .filter { sample_meta, new_refs -> [sample_meta, new_refs != ''] }
    //    .map { sample_meta, new_refs -> [sample_meta, new_refs.split('\t')] }
    //    .map { sample_meta, new_refs -> [sample_meta] }
    //    .unique() // sample_meta, ref_meta, ref_group_id

    //// Add new reference accessions to sample and reference metadata

    //sample_data = sample_data
    //    .transpose(by: 3)
    //    .combine(reference_data, by: 3)
    //    .map { ref_group_id, sample_meta, sample_paths, ncbi_accession, report_meta, ref_meta, ref_path, ref_ncbi_accession ->
    //         [sample_meta, sample_paths, ncbi_accession, report_meta, ref_meta]
    //    }
    //    .groupTuple(by: 0..3)

    //// Download reference files if an accession is provided
    //ref_ncbi_acc = reference_data
    //    .filter { ref_meta, path, ncbi_accession, ref_group_id ->
    //        ncbi_accession != null
    //    }
    //    .map { ref_meta, path, ncbi_accession, ref_group_id ->
    //        [ref_meta, ncbi_accession]
    //    }
    //    .unique()
    //DOWNLOAD_ASSEMBLIES ( ref_ncbi_acc )
    //reference_data = reference_data
    //    .join(DOWNLOAD_ASSEMBLIES.out.sequence, remainder: true)
    //    .map { ref_meta, path, ncbi_accession, ref_group_id, assembly ->
    //        [ref_meta, assembly ?: path]
    //    }


    emit:
    sample_data       // sample_meta, [sample_paths], report_meta, [ref_meta]
    reference_data    // ref_meta, path
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
    def sample_meta = [
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
    def report_meta = [id: row.report_group_ids]
    return [sample_meta, paths, ncbi_accession, ref_group_ids, report_meta]
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


