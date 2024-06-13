include { REFERENCE_INDEX        } from './reference_index'  // this is being called from subworkflows
include { ALIGN_READS            } from './align_reads'
include { CALL_VARIANTS          } from './call_variants'
include { IQTREE2 as IQTREE2_SNP } from '../../modules/local/iqtree2'
include { VCF_TO_TAB             } from '../../modules/local/vcf_to_tab'
include { VCF_TO_SNPALN          } from '../../modules/local/vcf_to_snpaln'
include { SEQKIT_SLIDING         } from '../../modules/nf-core/seqkit/sliding/main'
include { ASSIGN_MAPPING_REFERENCE                    } from '../../modules/local/assign_mapping_reference'

workflow VARIANT_ANALYSIS {

    take:
    sample_data
    ani_matrix // report_group_id, ani_matrix

    main:
    versions = Channel.empty()
    messages = Channel.empty()

    // Make file with smaple IDs and user-defined references or NA for each group
    samp_ref_pairs = sample_data
        .map{ [[id: it.sample_id], [id: it.report_group_id], it.ref_metas] }
        .transpose(by: 2)
        .map{ sample_id, report_id, ref_meta ->
            [sample_id, report_id, [id: ref_meta.ref_id], ref_meta.ref_path]
        }
        .tap{ references }
        .collectFile() { sample_id, report_id, ref_id, ref_path ->
            [ "${report_id.id}.csv", "${sample_id.id},${ref_id.id}\n" ]
        }
        .map {[[id: it.getSimpleName()], it]} // TODO this recreates the group_meta, but if other feilds besids "id" are added this will not preserve those

    // For each group, assign references for variant calling if not user-defined
    ASSIGN_MAPPING_REFERENCE (
        ani_matrix.join(samp_ref_pairs),
        params.ref_min_ani
    )
    mapping_ref_ids = ASSIGN_MAPPING_REFERENCE.out.samp_ref_pairs
        .splitText( elem: 1 )
        .map { [it[0], it[1].replace('\n', '')] } // remove newline that splitText adds
        .splitCsv( elem: 1 )
        .map { group_meta, sample_id, ref_id ->
            [[id: sample_id], group_meta, [id: ref_id]]
        }
        .join(references)

    //// Cutting up long reads
    //longreads = sample_data
    //    .filter { it.sequence_type != "illumina" }
    //SEQKIT_SLIDING ( longreads.map {  } )
    //chopped_reads = SEQKIT_SLIDING.out.fastx // meta, [chopped_reads]
    //    .join(longreads) // meta, [chopped_reads], [reads], ref_meta, reference, group_meta
    //    .map { it[0..1] + it[3..5] } // meta, [fastqs], ref_meta, reference, group_meta
    //sample_data = input
    //    .filter { it[0].reads_type == "illumina" }
    //    .mix(chopped_reads) // meta, [fastqs], ref_meta, reference, group_meta

    //// Remove any samples that do not have reference information
    //sample_data
    //    .branch {
    //        filtered: it[3] != null
    //        no_ref: it[3] == null
    //    }
    //    .set { parsed_sample_data }

    //// Report samples that do not have reference information
    //no_ref_warnings = parsed_sample_data.no_ref
    //    .map { [it[0], it[4], null, "VARIANT_ANALYSIS", "WARNING", "Sample is excluded from variant calling analysis because no reference genome is available."] } // meta, group_meta, ref_meta, workflow, level, message
    //messages = messages.mix(no_ref_warnings)

    //ch_ref = parsed_sample_data.filtered
    //    .map { it[2..3] }
    //    .groupTuple()
    //    .map { [it[0], it[1].sort()[0]] }
    //REFERENCE_INDEX ( parsed_sample_data.filtered.map { it[2..3] }.unique() )
    //versions = versions.mix(REFERENCE_INDEX.out.versions)

    //input_with_indexes = parsed_sample_data.filtered
    //    .map { [it[2], it[0], it[1], it[3], it[4]] } // [val(ref_meta), val(meta), [file(fastq)], file(reference), val(group_meta)]
    //    .combine(REFERENCE_INDEX.out.samtools_fai, by: 0)
    //    .combine(REFERENCE_INDEX.out.bwa_index, by: 0)
    //    .combine(REFERENCE_INDEX.out.picard_dict, by: 0)
    //    .map { [it[1], it[2], it[0], it[3], it[4], it[5], it[6], it[7]] } // [val(meta), [file(fastq)], val(ref_meta), file(reference), val(group_meta), fai, bwa, picard]

    //ALIGN_READS (
    //    input_with_indexes
    //        .map { it[0..3] + it[5..6] }
    //        .unique()
    //)
    //versions = versions.mix(ALIGN_READS.out.versions)

    //CALL_VARIANTS (
    //    input_with_indexes
    //        .map { [it[0], it[2], it[1]] + it[3..5] + [it[7]] } // [val(meta), val(ref_meta), [file(fastq)], file(reference), val(group_meta), fai, picard]
    //        .combine(ALIGN_READS.out.bam, by: 0..1)
    //        .combine(ALIGN_READS.out.bai, by: 0..1) // [val(meta), val(ref_meta), [file(fastq)], file(reference), val(group_meta), fai, picard, bam, bai]
    //        .map { [it[0], it[7], it[8], it[1]] + it[3..6] }
    //)
    //versions = versions.mix(CALL_VARIANTS.out.versions)

    //VCF_TO_TAB ( CALL_VARIANTS.out.vcf )
    //versions = versions.mix(VCF_TO_TAB.out.versions.toSortedList().map{it[0]})

    //VCF_TO_SNPALN ( VCF_TO_TAB.out.tab )
    //versions = versions.mix(VCF_TO_SNPALN.out.versions.toSortedList().map{it[0]})

    //// Dont make trees for groups with less than 3 samples
    //align_with_samp_meta = VCF_TO_SNPALN.out.fasta // val(ref+group_meta), fasta
    //    .combine(CALL_VARIANTS.out.samples, by:0) // val(ref+group_meta), fasta, [meta]
    //align_for_tree = align_with_samp_meta
    //    .filter { it[2].size() >= 3 }
    //    .map { it[0..1] } // val(ref+group_meta), fasta
    //too_few_samp_warnings = align_with_samp_meta // val(ref+group_meta), fasta, [meta]
    //    .filter { it[2].size() < 3 }
    //    .map { [null, it[0].group, it[0].ref, "VARIANT_ANALYSIS", "WARNING", "Not enough samples to build a SNP tree."] } // meta, group_meta, ref_meta, workflow, level, message
    //messages = messages.mix(too_few_samp_warnings)

    //IQTREE2_SNP ( align_for_tree, [] )
    //versions = versions.mix(IQTREE2_SNP.out.versions)

    //phylogeny = IQTREE2_SNP.out.phylogeny.map{ [it[0].group, it[0].ref, it[1]] } // group_meta, ref_meta, tree
    //snp_align = VCF_TO_SNPALN.out.fasta.map{ [it[0].group, it[0].ref, it[1]] }
    //vcf = CALL_VARIANTS.out.vcf.map{ [it[0].group, it[0].ref, it[1]] }

    //results = CALL_VARIANTS.out.vcf // ref+group_meta, vcf
    //    .combine(VCF_TO_SNPALN.out.fasta, by:0) // ref+group_meta, vcf, align
    //    .join(IQTREE2_SNP.out.phylogeny, remainder:true, by:0) // ref+group_meta, vcf, align, tree
    //    .map{ [it[0].group, it[0].ref] + it[1..3] } // group_meta, ref_meta, vcf, align, tree

    emit:
    //picard_dict  = REFERENCE_INDEX.out.picard_dict
    //samtools_fai = REFERENCE_INDEX.out.samtools_fai
    //samtools_gzi = REFERENCE_INDEX.out.samtools_gzi
    //bwa_index    = REFERENCE_INDEX.out.bwa_index
    //phylogeny    = phylogeny   // group_meta, ref_meta, tree
    //snp_align    = snp_align   // group_meta, ref_meta, fasta
    //vcf          = vcf         // group_meta, ref_meta, fasta
    //results      = results     // group_meta, ref_meta, vcf, align, tree
    versions     = versions // versions.yml
    messages     = messages    // meta, group_meta, ref_meta, workflow, level, message

}

