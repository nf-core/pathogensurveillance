include { REFERENCE_INDEX          } from '../reference_index'
include { ALIGN_READS              } from '../align_reads'
include { CALL_VARIANTS            } from '../call_variants'
include { IQTREE as IQTREE_SNP     } from '../../../modules/nf-core/iqtree'
include { VCF_TO_SNP_ALIGN         } from '../../../modules/local/vcf_to_snp_align'
include { SEQKIT_SLIDING           } from '../../../modules/nf-core/seqkit/sliding'
include { ASSIGN_MAPPING_REFERENCE } from '../../../modules/local/assign_mapping_reference'

workflow VARIANT_ANALYSIS {

    take:
    original_sample_data
    ani_matrix

    main:
    versions = channel.empty()
    messages = channel.empty()

    // Remove samples belonging to groups with only one sample
    grouped_sample_data = original_sample_data
        .map{ sample_meta -> [[id: sample_meta.report_group_ids], sample_meta] }
        .groupTuple(by: 0)
    sample_data = grouped_sample_data
        .filter{ groupEntry -> groupEntry[1].size() > 1 }
        .transpose(by: 1)
        .map{ groupEntry -> groupEntry[1] }
    messages = messages.mix(
        grouped_sample_data
            .filter{ groupEntry -> groupEntry[1].size() == 1 }
            .map{ report_meta, samp_metas ->
                [[id: samp_metas[0].sample_id], report_meta, null, "VARIANT_ANALYSIS", "WARNING", "Sample is excluded from variant calling analysis because it is the only sample in its report group."]
            }
    )

    // Make file with sample IDs and user-defined references or NA for each group
    samp_ref_pairs = sample_data
        .map{ sample_meta -> [sample_meta.sample_id, sample_meta.report_group_ids, sample_meta.ref_metas] }
        .transpose(by: 2)
        .map{ sample_id, report_group_id, ref_meta ->
            [sample_id, report_group_id, ref_meta.ref_id, ref_meta.ref_name, ref_meta.ref_description, ref_meta.ref_path, ref_meta.ref_primary_usage]
        }
        .unique()
        .tap{ references }
        .collectFile() { sample_id, report_group_id, ref_id, ref_name, ref_desc, ref_path, usage ->
            [ "${report_group_id}.tsv", "${sample_id}\t${ref_id}\t${ref_name}\t${ref_desc}\t${usage}\n" ]
        }
        .map {path -> [[id: path.getSimpleName()], path]}

    // For each group, assign references for variant calling if not user-defined
    ASSIGN_MAPPING_REFERENCE (
        ani_matrix.join(samp_ref_pairs),
        params.ref_min_ani
    )
    versions = versions.mix(ASSIGN_MAPPING_REFERENCE.out.versions)
    ref_paths = references
        .map {sample_id, report_group_id, ref_id, ref_name, ref_desc, ref_path, usage ->
            [[id: sample_id], [id:report_group_id], [id: ref_id], ref_path, usage]
        }
        .unique()
    sample_data_with_refs = ASSIGN_MAPPING_REFERENCE.out.samp_ref_pairs
        .splitText( elem: 1 )
        .map { row -> [row[0], row[1].replace('\n', '')] } // remove newline that splitText adds
        .splitCsv( elem: 1, sep: '\t' )
        .map { report_meta, tsv_contents ->
            [[id: tsv_contents[0]], report_meta, [id: tsv_contents[1]]]
        }
        .join(ref_paths, by: 0..2)
        .join(sample_data.map{ sample_meta -> [[id: sample_meta.sample_id], [id: sample_meta.report_group_ids], sample_meta.paths, sample_meta.sequence_type, sample_meta.ploidy] }, by: 0..1)
        .branch { entry -> // Remove any samples that do not have reference information
            filtered: entry[2] != null
            no_ref: entry[2] == null
        } // sample_meta, report_meta, ref_meta, ref_path, usage, read_paths, sequence_type

    // Remove samples belonging to groups with only one sample
    grouped_sample_data_with_refs = sample_data_with_refs.filtered
        .map{sample_meta, report_meta, ref_meta, ref_path, usage, read_paths, sequence_type, ploidy ->
            [report_meta, ref_meta, [sample_meta, report_meta, ref_meta, ref_path, usage, read_paths, sequence_type, ploidy]]
        }
        .groupTuple(by: [0,1])
    filtered_sample_data_with_refs = grouped_sample_data_with_refs
        .filter{ groupEntry -> groupEntry[2].size() >= 2 }
        .transpose(by: 2)
        .map{ groupEntry -> groupEntry[2] }
    messages = messages.mix(
        grouped_sample_data_with_refs
            .filter{ groupEntry -> groupEntry[2].size() < 2 }
            .map{ report_meta, ref_meta, data ->
                [data[0][0], report_meta, ref_meta, "VARIANT_ANALYSIS", "WARNING", "Sample is excluded from variant calling analysis because there are too few samples aligned to this reference to make a tree."]
            }
    )

    // Cutting up long reads
    longreads = filtered_sample_data_with_refs
        .filter { sample_meta, report_meta, ref_meta, ref_path, usage, read_paths, sequence_type, ploidy ->
            sequence_type == "nanopore" || sequence_type == "pacbio"
        }
    SEQKIT_SLIDING (
        longreads.map { sample_meta, report_meta, ref_meta, ref_path, usage, read_paths, sequence_type, ploidy ->
            [sample_meta, read_paths]
        }
        .unique()
    )
    chopped_reads = SEQKIT_SLIDING.out.fastx
        .combine(longreads, by: 0)
        .map { sample_meta, chopped_reads, report_meta, ref_meta, ref_path, usage, read_paths, sequence_type, ploidy ->
            [sample_meta, report_meta, ref_meta, ref_path, usage, chopped_reads, sequence_type, ploidy]
        }
    filtered_input = filtered_sample_data_with_refs
        .filter { sample_meta, report_meta, ref_meta, ref_path, usage, read_paths, sequence_type, ploidy ->
            sequence_type == "illumina" || sequence_type == "bgiseq"
        }
        .mix(chopped_reads) // meta, [fastqs], ref_meta, reference, report_meta


    // Report samples that do not have reference information
    no_ref_warnings = sample_data_with_refs.no_ref
        .map{ sample_meta, report_meta, ref_meta, ref_path, usage, read_paths, sequence_type, ploidy ->
            [sample_meta, report_meta, null, "VARIANT_ANALYSIS", "WARNING", "Sample is excluded from variant calling analysis because no reference genome is available."]
        }
    messages = messages.mix(no_ref_warnings)

    // Create indexes for each reference
    REFERENCE_INDEX (
        filtered_input
            .map{ sample_meta, report_meta, ref_meta, ref_path, usage, read_paths, sequence_type, ploidy ->
                [ref_meta, ref_path]
            }
            .unique()
    )
    versions = versions.mix(REFERENCE_INDEX.out.versions)

    input_with_indexes = filtered_input
        .map{ sample_meta, report_meta, ref_meta, ref_path, usage, read_paths, sequence_type, ploidy ->
            [ref_meta, sample_meta, read_paths, ref_path, report_meta]
        }
        .combine(REFERENCE_INDEX.out.samtools_fai, by: 0)
        .combine(REFERENCE_INDEX.out.bwa_index, by: 0)
        .combine(REFERENCE_INDEX.out.picard_dict, by: 0)
        .map{ ref_meta, sample_meta, read_paths, ref_path, report_meta, fai, bwa, picard ->
            [sample_meta, read_paths, ref_meta, ref_path, report_meta, fai, bwa, picard]
        }

    ALIGN_READS (
        input_with_indexes
            .map { row -> row[0..3] + row[5..6] }
            .unique()
    )
    versions = versions.mix(ALIGN_READS.out.versions)

    CALL_VARIANTS (
        input_with_indexes
            .map { row -> [row[0], row[2], row[1]] + row[3..5] + [row[7]] } // [val(meta), val(ref_meta), [file(fastq)], file(reference), val(report_meta), fai, picard]
            .combine(ALIGN_READS.out.bam, by: 0..1)
            .combine(ALIGN_READS.out.csi, by: 0..1) // [val(meta), val(ref_meta), [file(fastq)], file(reference), val(report_meta), fai, picard, bam, csi]
            .map { row -> [row[0], row[7], row[8], row[1]] + row[3..6] }
    )
    versions = versions.mix(CALL_VARIANTS.out.versions)

    sample_ploidy_data = filtered_sample_data_with_refs
        .collectFile(keepHeader: true, skip: 1) { sample_meta, report_meta, ref_meta, ref_path, usage, read_paths, sequence_type, ploidy ->
            [ "${report_meta.id}--${ref_meta.id}.tsv", "mapping_id\tploidy\n${ref_meta.id}--${sample_meta.id}\t${ploidy}\n" ]
        }
        .map { path -> [path.getSimpleName(), path] }
        .combine(CALL_VARIANTS.out.vcf.map { combined_meta, vcf -> [combined_meta.id, combined_meta]}, by: 0) // add on full combined metadata uing combined ID
        .map { combined_id, ploidy_data_path, combined_meta ->
            [combined_meta, ploidy_data_path]
        }

    VCF_TO_SNP_ALIGN (
        CALL_VARIANTS.out.vcf
            .combine(sample_ploidy_data, by: 0),
        params.max_variants
    )
    versions = versions.mix(VCF_TO_SNP_ALIGN.out.versions)
    removed_samps = VCF_TO_SNP_ALIGN.out.removed_sample_ids
        .splitText()
        .map { row -> [[id: row[1].replace('\n', '')], row[0].group, row[0].ref, "VARIANT_ANALYSIS", "WARNING", "Sample removed from SNP phylogeny due to too much missing data."] } // meta, group_meta, ref_meta, workflow, level, message

    // Dont make trees for groups with less than 3 samples
    align_with_samp_meta = VCF_TO_SNP_ALIGN.out.fasta
        .combine(VCF_TO_SNP_ALIGN.out.seq_count, by: 0)
        .branch { meta, fasta, seq_count_file ->
            enough: seq_count_file.text.trim().toInteger() >= 4
                return [meta, fasta]
            too_few: true
                return [meta, fasta]
        }
    too_few_samp_warnings = align_with_samp_meta.too_few
        .map { meta, fasta -> [null, meta.group, meta.ref, "VARIANT_ANALYSIS", "WARNING", "Not enough samples to build a SNP tree."] }
    messages = messages.mix(too_few_samp_warnings)

    phylogeny_input = align_with_samp_meta.enough
        .map { meta, alignment ->
            [meta, alignment, []]
        }
    IQTREE_SNP ( phylogeny_input, [], [], [], [], [], [], [], [], [], [], [], [] )
    phylogeny = IQTREE_SNP.out.phylogeny.map{ meta -> [meta[0].group, meta[0].ref, meta[1]] } // report_meta, ref_meta, tree
    snp_align = VCF_TO_SNP_ALIGN.out.fasta.map{ meta -> [meta[0].group, meta[0].ref, meta[1]] }
    vcf = CALL_VARIANTS.out.vcf.map{ meta -> [meta[0].group, meta[0].ref, meta[1]] }

    results = CALL_VARIANTS.out.vcf // ref+report_meta, vcf
        .combine(VCF_TO_SNP_ALIGN.out.fasta, by:0) // ref+report_meta, vcf, align
        .join(IQTREE_SNP.out.phylogeny, remainder:true, by:0) // ref+report_meta, vcf, align, tree
        .map{ meta -> [meta[0].group, meta[0].ref] + meta[1..3] } // report_meta, ref_meta, vcf, align, tree

    emit:
    picard_dict   = REFERENCE_INDEX.out.picard_dict
    samtools_fai  = REFERENCE_INDEX.out.samtools_fai
    samtools_gzi  = REFERENCE_INDEX.out.samtools_gzi
    bwa_index     = REFERENCE_INDEX.out.bwa_index
    phylogeny     = phylogeny   // report_meta, ref_meta, tree
    snp_align     = snp_align   // report_meta, ref_meta, fasta
    vcf           = vcf         // report_meta, ref_meta, fasta
    results       = results     // report_meta, ref_meta, vcf, align, tree
    mapping_ref   = ASSIGN_MAPPING_REFERENCE.out.samp_ref_pairs // report_meta, tsv
    versions      = versions // versions.yml
    messages      = messages    // meta, report_meta, ref_meta, workflow, level, message

}
