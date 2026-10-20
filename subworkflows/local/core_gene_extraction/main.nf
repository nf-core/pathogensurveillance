include { PIRATE                  } from '../../../modules/nf-core/pirate'
include { REFORMAT_PIRATE_RESULTS } from '../../../modules/local/reformat_pirate_results'
include { EXTRACT_PIRATE_GENES    } from '../../../modules/local/extract_pirate_genes'
include { CALCULATE_POCP          } from '../../../modules/local/calculate_pocp'
include { MAKE_GFF_WITH_FASTA     } from '../../../modules/local/make_gff_with_fasta'
include { BAKTA_BAKTA             } from '../../../modules/nf-core/bakta/bakta'
include { BAKTA_BAKTADBDOWNLOAD   } from '../../../modules/nf-core/bakta/baktadbdownload'
include { UNTAR                   } from '../../../modules/nf-core/untar'

workflow CORE_GENE_EXTRACTION {

    take:
    sample_data       // already domain-filtered and restricted to samples with assemblies
    sample_assemblies // sample_id, assembly_path
    selected_ref_data // ref_meta, report_meta, ref_path, ref_name, gff_path

    main:

    messages = channel.empty()

    // Split references by whether they already have a GFF
    selected_ref_data_has_gff = selected_ref_data.filter { ref_meta, report_meta, ref_path, ref_name, gff_path -> gff_path }
    selected_ref_data_no_gff = selected_ref_data.filter { ref_meta, report_meta, ref_path, ref_name, gff_path -> !gff_path }

    // Download the bakta database if needed
    //   Based on code from the bacass nf-core pipeline using the MIT license: https://github.com/nf-core/bacass
    if (params.bakta_db) {
        if (params.bakta_db.endsWith('.tar.gz')) {
            bakta_db_tar = channel.fromPath(params.bakta_db).map{ path -> [ [id: 'baktadb'], path] }
            UNTAR( bakta_db_tar )
            bakta_db = UNTAR.out.untar.map{ meta, db -> db }.first()
        } else {
            bakta_db = channel.fromPath(params.bakta_db).first()
        }
    } else if (params.download_bakta_db) {
        def bakta_cache = file("${params.data_dir}/bakta_db")
        if (bakta_cache.exists() && bakta_cache.listFiles()) {
            bakta_db = Channel.fromPath("${params.data_dir}/bakta_db/db-${params.bakta_db_type}").first()
        } else {
            BAKTA_BAKTADBDOWNLOAD()
            bakta_db = BAKTA_BAKTADBDOWNLOAD.out.db
        }
    }

    // Run bakta on samples and references without a gff already
    sample_assem_data = sample_data
        .map{ sample_meta -> [[id: sample_meta.sample_id], [id: sample_meta.report_group_ids]] }
        .combine(sample_assemblies, by: 0)
    ref_assem_data = selected_ref_data_no_gff
        .map { ref_meta, report_meta, ref_path, ref_name, gff_path ->
            [ref_meta, report_meta, ref_path]
        }
    all_assem_data = ref_assem_data
        .mix(sample_assem_data)
    BAKTA_BAKTA (
        all_assem_data
            .map { ref_meta, report_meta, assem_path ->
                [ref_meta, assem_path]
            }
            .unique(),
        bakta_db, // Bakta database
        [], // proteins (optional)
        [], // prodigal_tf (optional)
        [], // regions (optional)
        []  // hmms (optional)
    )

    // For references that have a gff already, combine the assembly with the gff the same way bakta outputs
    MAKE_GFF_WITH_FASTA (
        selected_ref_data_has_gff
            .map{ ref_meta, report_meta, ref_path, ref_name, ref_gff ->
                [ref_meta, ref_path, ref_gff]
            }
    )

    // group samples by report group
    bakta_gffs = all_assem_data
        .combine(BAKTA_BAKTA.out.gff, by: 0) // sample_or_ref_meta, report_meta, assem_path, gff_path
        .map { sample_or_ref_meta, report_meta, assem_path, gff_path ->
            [sample_or_ref_meta, report_meta, gff_path]
        }
    other_gffs = selected_ref_data_has_gff
        .combine(MAKE_GFF_WITH_FASTA.out.gff, by: 0)
        .map{ ref_meta, report_meta, ref_path, ref_name, ref_gff, ref_combined ->
            [ref_meta, report_meta, ref_combined]
        }
    all_gffs = bakta_gffs
        .mix(other_gffs)
        .map { sample_or_ref_meta, report_meta, gff ->
            [report_meta, gff]
        }
        .unique()
    PIRATE (
        all_gffs
            .groupTuple(by: 0, sort: 'hash')
            .map { meta, gffs -> // Attempt to fix intermittent "input file name collision" error
                def seen = [] as Set
                [meta, gffs.findAll { gff -> seen.add(gff.name) }]
            }
    )

    // Check that Pirate worked and report
    good_pirate_results = PIRATE.out.results
        .filter { result -> result[1].any{ path -> path.endsWith("PIRATE.gene_families.ordered.tsv") } }
    pirate_failed = PIRATE.out.results // val(group_meta), [result_files]
        .filter { result -> ! result[1].any{ path -> path.endsWith("PIRATE.gene_families.ordered.tsv") } }
        .map { result -> [null, result[0], null, "CORE_GENE_EXTRACTION", "WARNING", "Pirate failed to find a core genome, possibly because samples are very different or there are too few reads."] }
    messages = messages.mix(pirate_failed)

    REFORMAT_PIRATE_RESULTS ( good_pirate_results )

    // Calculate POCP from presence/absence matrix of genes
    CALCULATE_POCP (
        REFORMAT_PIRATE_RESULTS.out.gene_fam_pa
    )

    // Extract sequences of all genes into the unified gene_sequences format
    EXTRACT_PIRATE_GENES ( good_pirate_results )

    emit:
    gene_sequences = EXTRACT_PIRATE_GENES.out.gene_sequences
        .map { group_meta, gene_seqs -> [group_meta, gene_seqs, 'core', 'CORE_GENE_EXTRACTION'] }
    pocp           = CALCULATE_POCP.out.pocp
    messages       = messages

}
