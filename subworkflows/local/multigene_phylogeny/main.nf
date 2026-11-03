include { ASSIGN_CONTEXT_REFERENCES } from '../../../modules/local/assign_context_references'
include { CORE_GENE_EXTRACTION      } from '../core_gene_extraction'
include { BUSCO_GENE_EXTRACTION     } from '../busco_gene_extraction'
include { BUILD_GENE_TREES          } from '../build_gene_trees'

workflow MULTIGENE_PHYLOGENY {

    take:
    sample_data
    ani_matrix
    sample_assemblies

    main:

    // Drop samples without a successful assembly before splitting by domain
    successful_sample_data = sample_data
        .map{ sample_meta -> [[id: sample_meta.sample_id], sample_meta] }
        .combine(sample_assemblies, by: 0)
        .map{ sample_meta, sample_data_map, assembly_path -> sample_data_map }

    // Normalize domains once and give each report-group/domain pair its own selection key.
    phylogeny_sample_data = successful_sample_data
        .map { sample_meta ->
            def domain = sample_meta.domain == 'Eukaryota' ? 'Eukaryota' : 'Prokaryota'
            def selection_meta = [
                id: "${sample_meta.report_group_ids}--domain_${domain}",
                report_id: sample_meta.report_group_ids,
                domain: domain
            ]
            [sample_meta, selection_meta]
        }

    prokaryote_sample_data = phylogeny_sample_data
        .filter { sample_meta, selection_meta -> selection_meta.domain == 'Prokaryota' }
        .map { sample_meta, selection_meta -> sample_meta }
    eukaryote_sample_data = phylogeny_sample_data
        .filter { sample_meta, selection_meta -> selection_meta.domain == 'Eukaryota' }
        .map { sample_meta, selection_meta -> sample_meta }

    // Expand references once for both selection input and metadata matching.
    sample_reference_data = phylogeny_sample_data
        .flatMap { sample_meta, selection_meta ->
            sample_meta.ref_metas.collect { ref_meta ->
                [sample_meta.sample_id, selection_meta, ref_meta]
            }
        }

    // Build one sample/reference table for each report-group/domain selection.
    samp_ref_pairs = sample_reference_data
        .map { sample_id, selection_meta, ref_meta ->
            [selection_meta, "${sample_id}\t${ref_meta.ref_id}\t${ref_meta.ref_name}\t${ref_meta.ref_description}\t${ref_meta.ref_contextual_usage}"]
        }
        .unique()
        .groupTuple(by: 0, sort: 'hash')
        .map { selection_meta, lines ->
            [selection_meta, (lines.join('\n') + '\n').bytes.encodeBase64().toString()]
        }

    // Reuse each report group's ANI matrix for every domain represented in that group.
    ani_matrix_by_domain = phylogeny_sample_data
        .map { sample_meta, selection_meta -> [[id: sample_meta.report_group_ids], selection_meta] }
        .unique()
        .combine(ani_matrix, by: 0)
        .map { report_meta, selection_meta, matrix -> [selection_meta, matrix] }

    ASSIGN_CONTEXT_REFERENCES(
        ani_matrix_by_domain.combine(samp_ref_pairs, by: 0),
        params.n_ref_closest,
        params.n_ref_closest_named,
        params.n_ref_context
    )

    no_contextual_refs = ASSIGN_CONTEXT_REFERENCES.out.references
        .filter { selection_meta, tsv -> tsv.size() == 0 }
        .map { selection_meta, tsv ->
            [
                null,
                [id: selection_meta.report_id],
                null,
                'MULTIGENE_PHYLOGENY',
                'WARNING',
                "No contextual references were selected for the ${selection_meta.domain} domain in this group. Some outputs may be less informative."
            ]
        }

    // Match selected reference IDs to metadata using the same report/domain key.
    all_ref_data = sample_reference_data
        .map { sample_id, selection_meta, ref_meta ->
            [[id: ref_meta.ref_id], selection_meta, ref_meta.ref_path, ref_meta.ref_name, ref_meta.gff]
        }
        .groupTuple(by: [0, 1], sort: 'hash')
        .map { ref_meta, selection_meta, ref_paths, ref_names, gff_paths ->
            [ref_meta, selection_meta, ref_paths.unique()[0], ref_names.unique()[0], gff_paths.unique()[0]]
        }

    selected_ref_data = ASSIGN_CONTEXT_REFERENCES.out.references
        .splitText(elem: 1)
        .map { row -> [row[0], row[1].replace('\n', '')] }
        .splitCsv(elem: 1, sep: '\t')
        .map { selection_meta, tsv_contents -> [[id: tsv_contents[0].toString()], selection_meta] }
        .join(all_ref_data, by: [0, 1])
        .map { ref_meta, selection_meta, ref_path, ref_name, gff_path ->
            [ref_meta, [id: selection_meta.report_id], ref_path, ref_name, gff_path, selection_meta.domain]
        }

    prokaryote_selected_ref_data = selected_ref_data
        .filter { ref_meta, report_meta, ref_path, ref_name, gff_path, domain -> domain == 'Prokaryota' }
        .map { ref_meta, report_meta, ref_path, ref_name, gff_path, domain ->
            [ref_meta, report_meta, ref_path, ref_name, gff_path]
        }
    eukaryote_selected_ref_data = selected_ref_data
        .filter { ref_meta, report_meta, ref_path, ref_name, gff_path, domain -> domain == 'Eukaryota' }
        .map { ref_meta, report_meta, ref_path, ref_name, gff_path, domain ->
            [ref_meta, report_meta, ref_path, ref_name, gff_path]
        }

    selected_refs = ASSIGN_CONTEXT_REFERENCES.out.references
        .collectFile(storeDir: "${params.outdir}/reference_data/selected", sort: 'deep') { selection_meta, refs ->
            def ref_ids = refs.text.readLines().findAll { it.trim() }.unique().sort()
            def content = ref_ids ? ref_ids.join('\n') + '\n' : ''
            ["${selection_meta.report_id}_contextual_references.tsv", content]
        }
        .map { refs ->
            [[id: refs.name.replaceFirst(/_contextual_references\.tsv$/, '')], refs]
        }

    // Extract shared genes from prokaryote samples
    CORE_GENE_EXTRACTION(
        prokaryote_sample_data,
        sample_assemblies,
        prokaryote_selected_ref_data
    )

    // Extract shared genes from eukaryote samples
    BUSCO_GENE_EXTRACTION(
        eukaryote_sample_data,
        sample_assemblies,
        eukaryote_selected_ref_data
    )

    // Match each gene set to the contextual references selected for its domain.
    core_gene_sequences = CORE_GENE_EXTRACTION.out.gene_sequences
        .map { group_meta, gene_seqs, source, workflow_name ->
            [[id: "${group_meta.id}--domain_Prokaryota", report_id: group_meta.id, domain: 'Prokaryota'], gene_seqs, source, workflow_name]
        }
    busco_gene_sequences = BUSCO_GENE_EXTRACTION.out.gene_sequences
        .map { group_meta, gene_seqs, source, workflow_name ->
            [[id: "${group_meta.id}--domain_Eukaryota", report_id: group_meta.id, domain: 'Eukaryota'], gene_seqs, source, workflow_name]
        }
    all_gene_seqs = core_gene_sequences
        .mix(busco_gene_sequences)
        .combine(samp_ref_pairs, by: 0)
        .map { selection_meta, gene_seqs, source, workflow_name, samp_ref_pairs_raw ->
            [[id: selection_meta.report_id], gene_seqs, [source: source, samp_ref_pairs: samp_ref_pairs_raw, workflow_name: workflow_name]]
        }

    // Build trees from shared genes
    BUILD_GENE_TREES(all_gene_seqs)

    // Collect messages
    messages = channel.empty()
    messages = messages.mix(no_contextual_refs)
    messages = messages.mix(CORE_GENE_EXTRACTION.out.messages)
    messages = messages.mix(BUSCO_GENE_EXTRACTION.out.messages)
    messages = messages.mix(BUILD_GENE_TREES.out.messages)

    emit:
    selected_refs  = selected_refs
    pocp           = CORE_GENE_EXTRACTION.out.pocp
    phylogeny      = BUILD_GENE_TREES.out.phylogeny
    gene_lists     = BUILD_GENE_TREES.out.gene_lists
    messages       = messages

}
