//
// Check input samplesheets and convert to channels
//
include { BBMAP_SENDSKETCH       } from '../../../modules/nf-core/bbmap/sendsketch'
include { SAMPLESHEET_CHECK      } from '../../../modules/local/custom/samplesheet_check'
include { SRATOOLS_FASTERQDUMP   } from '../../../modules/local/sratools/fasterqdump'
include { INITIAL_CLASSIFICATION } from '../../../modules/local/custom/initial_classification'
include { DOWNLOAD_ASSEMBLIES    } from '../../../modules/local/custom/download_assemblies'
include { FIND_ASSEMBLIES        } from '../../../modules/local/custom/find_assemblies'
include { PICK_ASSEMBLIES        } from '../../../modules/local/custom/pick_assemblies'
include { SEQKIT_HEAD            } from '../../../modules/nf-core/seqkit/head'
include { COUNT_READS            } from '../../../modules/local/custom/count_reads'

workflow PREPARE_INPUT {
    take:
    sample_data_tsv
    reference_data_tsv

    main:

    // Initalize channel to accumulate information about software versions used
    versions = Channel.empty()
    messages = Channel.empty()

    // Parse input tables
    SAMPLESHEET_CHECK ( sample_data_tsv, reference_data_tsv )
    sample_data = SAMPLESHEET_CHECK.out.sample_data
        .splitCsv ( header:true, sep:'\t', quote:'"' )
        .map { create_sample_metadata_channel(it) }
        .filter { sample_meta -> sample_meta.enabled.toBoolean() }
    reference_data = SAMPLESHEET_CHECK.out.reference_data
        .splitCsv ( header:true, sep:'\t', quote:'"' )
        .map { create_reference_metadata_channel(it) }
        .filter { ref_meta -> ref_meta.ref_enabled.toBoolean() }
    messages = messages.mix (
        SAMPLESHEET_CHECK.out.message_data
            .splitCsv ( header:true, sep:'\t', quote:'"' )
            .map { [
                it.sample_id == '' ? null : [id: it.sample_id],
                it.report_group_id == '' ? null : [id: it.report_group_id],
                it.reference_id == '' ? null : [id: it.reference_id],
                "PREPARE_INPUT",
                it.message_type,
                it.description
            ] }
    )

    // Subest samples if max_samples is used
    if (params.max_samples) {
        sample_data = sample_data.take( params.max_samples )
    }

    // Add all of the reference metadata to the sample metadata
    sample_data_without_refs = sample_data
        .filter{ sample_meta -> sample_meta.ref_ids.size() == 0 }
        .map{ sample_meta -> [sample_meta, []] }
    sample_data = sample_data
        .map{ sample_meta -> [sample_meta.ref_ids, sample_meta] }
        .filter{ ref_ids, sample_meta -> ref_ids.size() > 0 }
        .transpose(by: 0)
        .combine(reference_data.map{ ref_meta -> [ref_meta.ref_id, ref_meta] }, by: 0)
        .map { ref_id, sample_meta, ref_meta -> [sample_meta, ref_meta] }
        .groupTuple(by: 0, sort: 'hash')
        .mix(sample_data_without_refs)

    // Download FASTQ files if an NCBI accession is provided
    sample_data_with_acc = sample_data
        .filter{ sample_meta, ref_metas -> sample_meta.ncbi_accession}
    sample_data_without_acc = sample_data
        .filter{ sample_meta, ref_metas -> ! sample_meta.ncbi_accession }
        .map { sample_meta, ref_metas -> [sample_meta, ref_metas]}
    ncbi_acc_sample_key = sample_data_with_acc
        .map{ sample_meta, ref_metas -> [[id: sample_meta.ncbi_accession], sample_meta, ref_metas] }
    ncbi_acc = ncbi_acc_sample_key
        .map { ncbi_acc_meta, sample_meta, ref_metas ->
            [ncbi_acc_meta, ncbi_acc_meta.id]
        }
        .unique()
    SRATOOLS_FASTERQDUMP ( ncbi_acc )
    versions = versions.mix(SRATOOLS_FASTERQDUMP.out.versions)
    sample_data = SRATOOLS_FASTERQDUMP.out.reads
        .combine(ncbi_acc_sample_key, by: 0)
        .map { ncbi_acc_meta, reads_path, sample_meta, ref_metas ->
            if (reads_path instanceof Collection) {
                if (reads_path.size() <= 2) {
                    sample_meta.paths = reads_path
                } else {  // if there are a mixture of single and paired end reads
                    sample_meta.paths = reads_path.findAll{it ==~ /^.+_[12]\..+$/}
                }
            } else {
                sample_meta.paths = [reads_path]
            }
            sample_meta.single_end = sample_meta.paths.size() == 1
            [sample_meta, ref_metas]
        }
        .mix(sample_data_without_acc)

    // Look up approximate taxonomic classifications
    BBMAP_SENDSKETCH (
        sample_data
            .map { sample_meta, ref_metas ->
                [[id: sample_meta.sample_id], sample_meta.paths]
            }
            .unique()
    )
    versions = versions.mix(BBMAP_SENDSKETCH.out.versions)

    // Parse results of sendsketch to get list of taxa to download references for
    INITIAL_CLASSIFICATION ( BBMAP_SENDSKETCH.out.hits )
    versions = versions.mix(INITIAL_CLASSIFICATION.out.versions)

    // Add estimated depth and kingdom to sample metadata
    sendsketch_depth = BBMAP_SENDSKETCH.out.hits
        .splitText(limit: 2, by: 2)
        .map { sample_meta, header ->
            def match = header =~ /Depth: ([0-9.]+)/
            [sample_meta, match[0][1]]
        }
    sample_data = sample_data
        .map{ sample_meta, ref_metas ->
            [[id: sample_meta.sample_id], sample_meta, ref_metas]
        }
        .combine(sendsketch_depth, by: 0)
        .combine(INITIAL_CLASSIFICATION.out.kingdom, by: 0)
        .map{ sample_id, sample_meta, ref_metas, depth, kingdom ->
            sample_meta.sendsketch_depth = depth
            sample_meta.kingdom = kingdom
            [sample_meta, ref_metas]
        }

    // Get list of families for all samples without exclusive references defined by the user
    def no_auto_contextual_refs = params.n_ref_closest == 0 && params.n_ref_closest_named == 0 && params.n_ref_context == 0
    all_families = sample_data
        .filter { sample_meta, ref_metas -> // Dont lookup assembly metadata if no references are to be dowloaded automatically (besides user-defined references)
            ! (params.n_ref_genera == 0 && params.n_ref_species == 0 && params.n_ref_strains == 0)
        }
        .filter { sample_meta, ref_metas -> // Dont lookup assembly metadata for samples that the user has defined exclusive references for
            ! (ref_metas.collect{it.ref_primary_usage}.contains('exclusive') &&
                (ref_metas.collect{it.ref_contextual_usage}.contains('exclusive') || no_auto_contextual_refs))
        }
        .map { sample_meta, ref_metas ->
            [[id: sample_meta.sample_id], sample_meta]
        }
        .combine(INITIAL_CLASSIFICATION.out.families, by: 0)
        .map{ sample_id, sample_meta, families -> [families] }
        .splitText()
        .flatten()
        .map { families -> families.replace('\n', '') }
        .unique()

    // Download RefSeq metadata for all assemblies for every family found by the initial identification
    FIND_ASSEMBLIES (
        all_families
    )
    versions = versions.mix(FIND_ASSEMBLIES.out.versions)

    // Add placeholders for NCBI reference metadata if none was looked up
    ncbi_ref_meta = INITIAL_CLASSIFICATION.out.families
        .splitText(elem: 1)
        .map { sample_meta, families ->
            [families.replace('\n', '')]
        }
        .unique()
        .join(FIND_ASSEMBLIES.out.stats.ifEmpty([null, null]), remainder: true)
        .filter { it != [null, null] }

    // Choose reference sequences to provide context for each sample
    family_stats_per_sample = INITIAL_CLASSIFICATION.out.families
        .splitText(elem: 1)
        .map { sample_id, families -> [families.replace('\n', ''), sample_id] }
        .combine(ncbi_ref_meta, by: 0)
        .map { family, sample_id, stats_path -> [sample_id, stats_path] }
        .unique()
        .groupTuple(by: 0, sort: 'hash')
        .map { sample_id, family_stats ->
            [sample_id, family_stats.findAll{it != null}]
        }
    taxon_and_ref_data = sample_data
        .map { sample_meta, ref_metas -> [[id: sample_meta.sample_id]] }
        .unique()
        .join(INITIAL_CLASSIFICATION.out.families)
        .join(INITIAL_CLASSIFICATION.out.genera)
        .join(INITIAL_CLASSIFICATION.out.species)
        .join(family_stats_per_sample)
        .filter { sample_meta, families, genera, species, family_stats ->
            family_stats.size() > 0
        }
    PICK_ASSEMBLIES (
        taxon_and_ref_data,
        params.n_ref_strains,
        params.n_ref_species,
        params.n_ref_genera,
        params.only_latin_binomial_refs
    )

    // Add placeholders for PICK_ASSEMBLIES output if not run
    picked_assemblies_stat_files = sample_data
        .map { sample_meta, ref_metas -> [[id: sample_meta.sample_id]] }
        .unique()
        .join(PICK_ASSEMBLIES.out.stats.ifEmpty([null, null]), remainder: true)
        .filter { it != [null, null] } // above join adds [null, null] if channel is empty
    picked_assemblies_refs = PICK_ASSEMBLIES.out.stats // pick_assemblies_out has a list of ref metdata for each sample
        .splitCsv(header:true, sep:'\t', quote:'"', elem: 1)
        .map{ sample_id, ref_meta ->
            [sample_id, ref_meta.collectEntries{ key, value -> [(key): value ?: null] }]
        }
        .unique()
        .groupTuple(by: 0, sort: 'hash')
    new_reference_data = sample_data
        .map { sample_meta, ref_metas -> [[id: sample_meta.sample_id]] }
        .unique()
        .join(picked_assemblies_refs.ifEmpty([null, null]), remainder: true)
        .filter { it != [null, null] } // above join adds [null, null] if channel is empty
        .map { sample_meta, ref_metas ->
            [sample_meta, ref_metas == null ? [] : ref_metas]
        }

    // Add new refernces to sample metadata
    sample_data = sample_data
        .map { sample_meta, ref_metas ->
            [[id: sample_meta.sample_id], sample_meta, ref_metas]
        }
        .combine(new_reference_data, by: 0)
        .map { sample_id, sample_meta, ref_metas_user, ref_metas_to_download ->
            [sample_id, ref_metas_user + (ref_metas_to_download ?: [])]
        }
        .combine(sample_data.map{ sample_meta, ref_metas -> [[id: sample_meta.sample_id], sample_meta] }, by: 0)
        .map { sample_id, ref_metas, sample_meta ->
            [sample_meta, ref_metas]
        }

    reference_data = new_reference_data
        .transpose(by: 1)
        .map{ sample_id, ref_meta -> ref_meta }
        .mix(reference_data)

    // Warn of samples for which no reference information could be found
    no_assemblies_found = PICK_ASSEMBLIES.out.line_count
        .filter { sample_id, line_count ->
            line_count == "1"
        }
        .combine(sample_data.map { sample_meta, ref_metas -> [[id: sample_meta.sample_id], sample_meta, ref_metas] }, by: 0)
        .map { sample_id, line_count, sample_meta, ref_metas ->
            [sample_meta, ref_metas]
        }
    no_ref_warnings = no_assemblies_found
        .map{ sample_meta, ref_metas ->
            [sample_meta, [id: sample_meta.report_group_ids], null, "PREPARE_INPUT", "WARNING", "Could not find any references to download."]
        }
    messages = messages.mix(no_ref_warnings)

    // Download reference files if an accession is provided
    ref_ncbi_acc = reference_data
        .filter{ ref_meta -> ref_meta.ref_ncbi_accession }
        .tap{ ref_data_with_ncbi_acc }
        .map { ref_meta ->
            [[id: ref_meta.ref_ncbi_accession], ref_meta.ref_ncbi_accession]
        }
        .unique()
    DOWNLOAD_ASSEMBLIES ( ref_ncbi_acc )
    versions = versions.mix(DOWNLOAD_ASSEMBLIES.out.versions)
    local_references = sample_data
        .transpose(by: 1)
        .filter{ sample_meta, ref_meta ->
            ref_meta.ref_path
        }
    sample_data = sample_data
        .transpose(by: 1)
        .map{ sample_meta, ref_meta ->
            [[id: ref_meta.ref_ncbi_accession], sample_meta, ref_meta ]
        }
        .combine(DOWNLOAD_ASSEMBLIES.out.sequence, by: 0)
        .combine(DOWNLOAD_ASSEMBLIES.out.gff, by: 0)
        .map { ncbi_acc_meta, sample_meta, ref_meta, ref_path, gff_path ->
            ref_meta.ref_path = ref_path
            ref_meta.gff = gff_path
            [sample_meta, ref_meta]
        }
        .mix(local_references)
        .unique()
        .groupTuple(by: 0, sort: 'hash')

    // Add reference metadata list to the sample metadata
    sample_data = sample_data
        .map{ sample_meta, ref_metas ->
            sample_meta.ref_metas = ref_metas
            sample_meta
        }

    // Ensure that items that can be single or multiple are always formatted as lists
    sample_data = sample_data
        .map{ sample_meta ->
            sample_meta.paths = sample_meta.paths instanceof Collection ? sample_meta.paths : [sample_meta.paths]
            sample_meta
        }

    // Count the number of reads and basepairs to decide whether not to subset_reads
    samples_to_subset = sample_data
        .map { [[id: it.sample_id], it.paths, it.sendsketch_depth] }
        .unique()
        .filter { sample_id, fastq_paths, depth ->
            depth.toFloat() > params.max_depth.toFloat()
        }
    samples_to_not_subset = sample_data
        .filter {
            it.sendsketch_depth.toFloat() <= params.max_depth.toFloat()
        }
    COUNT_READS (
        samples_to_subset
            .map { sample_meta, fastq_paths, depth ->
                [sample_meta, fastq_paths]
            }
            .unique(),
    )

    // Subset sample reads to increase speed of following steps
    SEQKIT_HEAD (
        samples_to_subset
            .combine(COUNT_READS.out.read_count, by: 0)
            .map { sample_meta, fastq_paths, depth, read_count ->
                [sample_meta, fastq_paths, Math.ceil((params.max_depth.toFloat() / depth.toFloat()) * read_count.toFloat()).toInteger() ]
            }
    )
    versions = versions.mix(SEQKIT_HEAD.out.versions)
    sample_data = sample_data
        .map { sample_meta ->
            [[id: sample_meta.sample_id], sample_meta]
        }
        .combine(SEQKIT_HEAD.out.subset, by: 0)
        .map { sample_id, sample_meta, subset_reads ->
            sample_meta.paths = subset_reads
            sample_meta
        }
        .mix(samples_to_not_subset)

    emit:
    sample_data
    sendsketch = BBMAP_SENDSKETCH.out.hits
    family_stats = ncbi_ref_meta
    selected_ref_meta = picked_assemblies_stat_files
    family_stats_per_sample = family_stats_per_sample
    versions = SAMPLESHEET_CHECK.out.versions
    messages = messages    // meta, group_meta, ref_meta, workflow, level, message
}

def create_sample_metadata_channel(LinkedHashMap sample_meta) {
    if (sample_meta.path != "" && !file(sample_meta.path).exists()) {
        exit 1, "ERROR: Please check the sample metadata TSV/CSV. The file specified by 'path` does not exist.\n${sample_meta.path}"
    }
    if (sample_meta.path_2 != "" && !file(sample_meta.path_2).exists()) {
        exit 1, "ERROR: Please check the sample metadata TSV/CSV. The file specified by 'path_2' does not exist.\n${sample_meta.path_2}"
    }
    sample_meta = sample_meta.collectEntries { key, value -> [(key): value ?: null] }
    sample_meta.ref_ids = sample_meta.ref_ids ? sample_meta.ref_ids.split(";") as ArrayList : []
    sample_meta.single_end = ! sample_meta.path_2
    def paths = null
    if (sample_meta.path) {
        paths = [file(sample_meta.path)]
        if (sample_meta.path_2) {
            paths.add(file(sample_meta.path_2))
        }
    }
    sample_meta.paths = paths
    return sample_meta
}

def create_reference_metadata_channel(LinkedHashMap ref_meta) {
    if (ref_meta.ref_path != "" && !file(ref_meta.ref_path).exists()) {
        exit 1, "ERROR: Please check the reference metadata TSV/CSV. The file specified by 'ref_path` does not exist.\n${ref_meta.ref_path}"
    }
    ref_meta = ref_meta.collectEntries { key, value -> [(key): value ?: null] }
    ref_meta.ref_path = ref_meta.ref_path ? file(ref_meta.ref_path) : null
    return ref_meta
}


