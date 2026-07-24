//
// Check input samplesheets and convert to channels
//
include { BBMAP_SENDSKETCH       } from '../../../modules/nf-core/bbmap/sendsketch'
include { SAMPLESHEET_CHECK      } from '../../../modules/local/samplesheet_check'
include { SRATOOLS_FASTERQDUMP   } from '../../../modules/nf-core/sratools/fasterqdump'
include { INITIAL_CLASSIFICATION } from '../../../modules/local/initial_classification'
include { DOWNLOAD_ASSEMBLIES    } from '../../../modules/local/download_assemblies'
include { FIND_ASSEMBLIES        } from '../../../modules/local/find_assemblies'
include { PARSE_ASSEMBLIES       } from '../../../modules/local/parse_assemblies'
include { PICK_ASSEMBLIES        } from '../../../modules/local/pick_assemblies'
include { SEQKIT_HEAD            } from '../../../modules/nf-core/seqkit/head'
include { SEQKIT_STATS           } from '../../../modules/nf-core/seqkit/stats'

workflow PREPARE_INPUT {
    take:
    sample_data_tsv
    reference_data_tsv

    main:

    // Initalize channel to accumulate information about software versions used
    versions = channel.empty()
    messages = channel.empty()

    // Parse input tables
    SAMPLESHEET_CHECK ( sample_data_tsv, reference_data_tsv, params.max_samples )
    versions = versions.mix(SAMPLESHEET_CHECK.out.versions)
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

    // Check for cached reads and split into cached / to-download
    ncbi_acc_cached = ncbi_acc.map { ncbi_acc_meta, sra ->
        def cache_dir = file("${params.data_dir}/reads")
        def cached = cache_dir.exists() ? cache_dir.list().findAll { filename -> filename.startsWith(ncbi_acc_meta.id) && filename.endsWith('.fastq.gz') } : []
        [ncbi_acc_meta, sra, cached]
    }
    ncbi_acc_with_cache = ncbi_acc_cached.filter { ncbi_acc_meta, sra, cached -> cached.size() > 0 }
    ncbi_acc_to_download = ncbi_acc_cached.filter { ncbi_acc_meta, sra, cached -> cached.size() == 0 }

    // Emit cached reads
    cached_reads = ncbi_acc_with_cache.map { ncbi_acc_meta, sra, cached ->
        def paths = cached.collect { file("${params.data_dir}/reads/${it}") }
        [ncbi_acc_meta, paths]
    }

    // Download missing reads
    SRATOOLS_FASTERQDUMP ( ncbi_acc_to_download.map { ncbi_acc_meta, sra, cached ->  [ncbi_acc_meta, sra] }, [], [] )
    downloaded_reads = SRATOOLS_FASTERQDUMP.out.reads

    // Combine cached and downloaded reads
    all_reads = cached_reads.mix(downloaded_reads)

    sample_data = all_reads
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

    // Extract the depth estimate from the sendsketch result and remove samples that dont have it
    sendsketch_depth = BBMAP_SENDSKETCH.out.hits
        .splitText(limit: 2, by: 2)
        .filter{sample_meta, header -> header =~ /Depth:/}
        .map { sample_meta, header ->
            def match = header =~ /Depth: ([0-9.]+)/
            [sample_meta, match[0][1]]
        }
    sendsketch_out = BBMAP_SENDSKETCH.out.hits
        .combine(sendsketch_depth, by: 0) // acts as a filter
        .map{ sample_meta, hit_data, depth -> [sample_meta, hit_data]}
    no_sketch_warnings =  BBMAP_SENDSKETCH.out.hits
        .splitText(limit: 2, by: 2)
        .filter{sample_meta, header -> header !=~ /Depth:/}
        .map{ sample_meta, empty ->
            [sample_meta, [id: sample_meta.report_group_ids], null, "PREPARE_INPUT", "WARNING", "Not enough data to make an initial classification. Check size of input data."]
        }
    messages = messages.mix(no_sketch_warnings)

    // Parse results of sendsketch to get list of taxa to download references for
    INITIAL_CLASSIFICATION ( sendsketch_out )
    versions = versions.mix(INITIAL_CLASSIFICATION.out.versions)

    // Add estimated depth and domain to sample metadata
    sample_data = sample_data
        .map{ sample_meta, ref_metas ->
            [[id: sample_meta.sample_id], sample_meta, ref_metas]
        }
        .combine(sendsketch_depth, by: 0)
        .combine(INITIAL_CLASSIFICATION.out.domain, by: 0)
        .map{ sample_id, sample_meta, ref_metas, depth, domain ->
            sample_meta.sendsketch_depth = depth
            sample_meta.domain = domain
            [sample_meta, ref_metas]
        }

    // Get list of families for all samples without exclusive references defined by the user
    family_taxon_ids = INITIAL_CLASSIFICATION.out.taxa_found
        .splitCsv(elem: 1, header: true, sep: '\t')
        .filter{ sample_meta, taxon_data ->
            taxon_data.rank == 'family'
        }
        .map{ sample_meta, taxon_data ->
            [sample_meta, taxon_data.taxon_id]
        }
    def no_auto_contextual_refs = params.n_ref_closest == 0 && params.n_ref_closest_named == 0 && params.n_ref_context == 0
    family_ids_to_download = sample_data
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
        .combine(family_taxon_ids, by: 0)
        .map{ sample_id, sample_meta, family_ids -> [family_ids] }
        .flatten()
        .unique()

    // Check for cached parsed assembly metadata
    family_ids_cached = family_ids_to_download.map { taxon ->
        def cache_file = file("${params.data_dir}/assembly_metadata/${taxon}.tsv")
        [taxon, cache_file.exists() ? cache_file : null]
    }
    family_ids_with_cache = family_ids_cached.filter { taxon, cached -> cached != null }
    family_ids_to_process = family_ids_cached.filter { taxon, cached -> cached == null }

    // Emit cached stats directly
    cached_stats = family_ids_with_cache.map { taxon, cached ->
        [taxon, cached]
    }

    // Download RefSeq metadata for all assemblies for every family found by the initial identification
    FIND_ASSEMBLIES (
        family_ids_to_process.map { taxon, cached -> taxon }
    )
    versions = versions.mix(FIND_ASSEMBLIES.out.versions)

    // Parse assembly metadata to TSVs to save time when multiple samples use the same data
    PARSE_ASSEMBLIES (
        FIND_ASSEMBLIES.out.stats
    )
    versions = versions.mix(PARSE_ASSEMBLIES.out.versions)
    processed_stats = PARSE_ASSEMBLIES.out.stats

    // Combine cached and processed stats
    all_stats = cached_stats.mix(processed_stats)

    // Add placeholders for NCBI reference metadata if none was looked up
    ncbi_ref_meta = family_taxon_ids
        .map { sample_meta, family_ids ->
            [family_ids]
        }
        .unique()
        .join(all_stats.ifEmpty([null, null]), remainder: true)
        .filter { item -> item != [null, null] }

    // Choose reference sequences to provide context for each sample
    family_stats_per_sample = family_taxon_ids
        .map { sample_id, family_id -> [family_id, sample_id] }
        .combine(ncbi_ref_meta, by: 0)
        .map { family_id, sample_id, stats_path -> [sample_id, stats_path] }
        .unique()
        .groupTuple(by: 0, sort: 'hash')
        .map { sample_id, family_stats ->
            [sample_id, family_stats.findAll{it != null}]
        }

    // Build a stable comma-delimited string of excluded accessions per sample
    excluded_accessions_per_sample = sample_data
        .map { sample_meta, ref_metas ->
            def excluded = ref_metas
                .findAll { ref_meta -> ref_meta.ref_primary_usage == 'excluded' && ref_meta.ref_contextual_usage == 'excluded' && ref_meta.ref_ncbi_accession }
                .collect { ref_meta -> ref_meta.ref_ncbi_accession }
                .sort()
                .join(',')
            [[id: sample_meta.sample_id], excluded ?: 'NONE']
        }

    taxon_and_ref_data = sample_data
        .map { sample_meta, ref_metas -> [[id: sample_meta.sample_id]] }
        .unique()
        .join(INITIAL_CLASSIFICATION.out.taxa_found)
        .join(family_stats_per_sample)
        .filter { sample_meta, taxa_found, family_stats ->
            family_stats.size() > 0
        }
        .join(excluded_accessions_per_sample)
        .map { sample_meta, taxa_found, family_stats, excluded_accessions ->
            [sample_meta, taxa_found, family_stats, excluded_accessions]
        }
    PICK_ASSEMBLIES (
        taxon_and_ref_data,
        params.n_ref_strains,
        params.n_ref_species,
        params.n_ref_genera,
        params.only_latin_binomial_refs
    )
    versions = versions.mix(PICK_ASSEMBLIES.out.versions)

    // Add placeholders for PICK_ASSEMBLIES output if not run
    picked_assemblies_stat_files = sample_data
        .map { sample_meta, ref_metas -> [[id: sample_meta.sample_id]] }
        .unique()
        .join(PICK_ASSEMBLIES.out.metadata.ifEmpty([null, null]), remainder: true)
        .filter { item -> item != [null, null] } // above join adds [null, null] if channel is empty
    picked_assemblies_refs = PICK_ASSEMBLIES.out.formatted // pick_assemblies_out has a list of ref metdata for each sample
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
        .filter { item -> item != [null, null] } // above join adds [null, null] if channel is empty
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
        .filter{ ref_meta -> ref_meta.ref_ncbi_accession && !(ref_meta.ref_primary_usage == 'excluded' && ref_meta.ref_contextual_usage == 'excluded') }
        .tap{ ref_data_with_ncbi_acc }
        .map { ref_meta ->
            [[id: ref_meta.ref_ncbi_accession], ref_meta.ref_ncbi_accession]
        }
        .unique()

    // Check for cached assemblies and split into cached / to-download
    ref_ncbi_acc_cached = ref_ncbi_acc.map { ref_meta, acc ->
        def fasta = file("${params.data_dir}/assemblies/${ref_meta.id}.fasta.gz")
        [ref_meta, acc, fasta.exists()]
    }
    ref_ncbi_acc_with_cache = ref_ncbi_acc_cached.filter { ref_meta, acc, has_fasta -> has_fasta }
    ref_ncbi_acc_to_download = ref_ncbi_acc_cached.filter { ref_meta, acc, has_fasta -> !has_fasta }

    // Emit cached assemblies
    cached_seq_and_gff = ref_ncbi_acc_with_cache.map { ref_meta, acc, has_fasta ->
        def fasta = file("${params.data_dir}/assemblies/${ref_meta.id}.fasta.gz")
        def gff = file("${params.data_dir}/assemblies/${ref_meta.id}.gff.gz")
        [ref_meta, fasta, gff.exists() ? gff : null]
    }

    // Download missing assemblies
    DOWNLOAD_ASSEMBLIES ( ref_ncbi_acc_to_download.map { ref_meta, acc, cached -> [ref_meta, acc] } )
    versions = versions.mix(DOWNLOAD_ASSEMBLIES.out.versions)
    downloaded_seq_and_gff = DOWNLOAD_ASSEMBLIES.out.sequence
        .join(DOWNLOAD_ASSEMBLIES.out.gff, by: 0, remainder: true)

    // Combine cached and downloaded assemblies
    all_seq_and_gff = cached_seq_and_gff.mix(downloaded_seq_and_gff)

    // Separate excluded references to preserve their metadata without downloading/processing
    excluded_refs = sample_data
        .transpose(by: 1)
        .filter{ sample_meta, ref_meta -> ref_meta.ref_primary_usage == 'excluded' && ref_meta.ref_contextual_usage == 'excluded' }
    local_references = sample_data
        .transpose(by: 1)
        .filter{ sample_meta, ref_meta ->
            ref_meta.ref_path && !(ref_meta.ref_primary_usage == 'excluded' && ref_meta.ref_contextual_usage == 'excluded')
        }
    sample_data = sample_data
        .transpose(by: 1)
        .filter{ sample_meta, ref_meta -> !(ref_meta.ref_primary_usage == 'excluded' && ref_meta.ref_contextual_usage == 'excluded') }
        .map{ sample_meta, ref_meta ->
            [[id: ref_meta.ref_ncbi_accession], sample_meta, ref_meta ]
        }
        .combine(all_seq_and_gff, by: 0)
        .map { ncbi_acc_meta, sample_meta, ref_meta, ref_path, gff_path ->
            ref_meta.ref_path = ref_path
            ref_meta.gff = gff_path
            [sample_meta, ref_meta]
        }
        .mix(local_references)
        .mix(excluded_refs)
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

    // Count the number of reads and basepairs to decide whether not to subset reads
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
    SEQKIT_STATS (
        samples_to_subset
            .map { sample_meta, fastq_paths, depth ->
                [sample_meta, fastq_paths]
            }
            .unique(),
    )
    read_count = SEQKIT_STATS.out.stats
        .splitCsv ( header:true, sep:'\t', limit: 1 )
        .map { sample_meta, stats ->
            [sample_meta, stats.sum_len]
        }

    // Subset sample reads to increase speed of following steps
    SEQKIT_HEAD (
        samples_to_subset
            .combine(read_count, by: 0)
            .map { sample_meta, fastq_paths, depth, my_read_count ->
                [sample_meta, fastq_paths, Math.ceil((params.max_depth.toFloat() / depth.toFloat()) * my_read_count.toFloat()).toInteger() ]
            }
    )
    versions = versions.mix(SEQKIT_HEAD.out.versions)
    sample_data = sample_data
        .map { sample_meta ->
            [[id: sample_meta.sample_id], sample_meta]
        }
        .combine(SEQKIT_HEAD.out.subset, by: 0)
        .map { sample_id, sample_meta, subset_reads ->
            sample_meta.paths = subset_reads instanceof Collection ? subset_reads : [subset_reads]
            sample_meta
        }
        .mix(samples_to_not_subset)

    emit:
    sample_data
    sendsketch = BBMAP_SENDSKETCH.out.hits
    taxa_found = INITIAL_CLASSIFICATION.out.taxa_found
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
