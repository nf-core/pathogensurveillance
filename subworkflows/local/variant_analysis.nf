include { REFERENCE_INDEX        } from './reference_index'  // this is being called from subworkflows
include { ALIGN_READS            } from './align_reads'
include { CALL_VARIANTS          } from './call_variants'
include { IQTREE2 as IQTREE2_SNP } from '../../modules/local/iqtree2'        
include { VCF_TO_TAB             } from '../../modules/local/vcf_to_tab'                   
include { VCF_TO_SNPALN          } from '../../modules/local/vcf_to_snpaln'       

workflow VARIANT_ANALYSIS {

    take:
    input // [val(meta), [file(fastq)], val(ref_meta), file(reference), val(group_meta)]
    ch_samplesheet // channel: path

    main:
    ch_versions = Channel.empty()
    messages = Channel.empty()
    
    // Remove any samples that do not have reference information
    input
        .branch {
            filtered: it[3] != null
            no_ref: it[3] == null
        }
        .set { parsed_input }
         
    
//    input_filtered = input
//        .filter { it[3] != null } // meta, [fastq], ref_meta, reference, group_meta
        
    // Report samples that do not have reference information
    no_ref_warnings = parsed_input.no_ref
        .map { [it[0], it[4], null, "VARIANT_ANALYSIS", "WARNING", "Sample is excluded from variant calling analysis because no reference genome is available."] } // meta, group_meta, ref_meta, workflow, level, message
    messages = messages.mix(no_ref_warnings)
    
    ch_ref = parsed_input.filtered
        .map { it[2..3] }
        .groupTuple()
        .map { [it[0], it[1].sort()[0]] }
    REFERENCE_INDEX ( parsed_input.filtered.map { it[2..3] }.unique() )
    ch_versions = ch_versions.mix(REFERENCE_INDEX.out.versions) 

    input_with_indexes = parsed_input.filtered
        .map { [it[2], it[0], it[1], it[3], it[4]] } // [val(ref_meta), val(meta), [file(fastq)], file(reference), val(group_meta)]
        .combine(REFERENCE_INDEX.out.samtools_fai, by: 0)              
        .combine(REFERENCE_INDEX.out.bwa_index, by: 0)
        .combine(REFERENCE_INDEX.out.picard_dict, by: 0)
        .map { [it[1], it[2], it[0], it[3], it[4], it[5], it[6], it[7]] } // [val(meta), [file(fastq)], val(ref_meta), file(reference), val(group_meta), fai, bwa, picard]     

    ALIGN_READS (
        input_with_indexes
            .map { it[0..3] + it[5..6] }
            .unique()
    )
    ch_versions = ch_versions.mix(ALIGN_READS.out.versions)       

    CALL_VARIANTS (
        input_with_indexes
            .map { [it[0], it[2], it[1]] + it[3..5] + [it[7]] } // [val(meta), val(ref_meta), [file(fastq)], file(reference), val(group_meta), fai, picard]
            .combine(ALIGN_READS.out.bam, by: 0..1)
            .combine(ALIGN_READS.out.bai, by: 0..1) // [val(meta), val(ref_meta), [file(fastq)], file(reference), val(group_meta), fai, picard, bam, bai]
            .map { [it[0], it[7], it[8], it[1]] + it[3..6] }
    )
    ch_versions = ch_versions.mix(CALL_VARIANTS.out.versions)       
   
    VCF_TO_TAB ( CALL_VARIANTS.out.vcf )                                                       
    ch_versions = ch_versions.mix(VCF_TO_TAB.out.versions.toSortedList().map{it[0]})
                                                                                
    VCF_TO_SNPALN ( VCF_TO_TAB.out.tab )                                            
    ch_versions = ch_versions.mix(VCF_TO_SNPALN.out.versions.toSortedList().map{it[0]})
    
    // Dont make trees for groups with less than 3 samples
    align_with_samp_meta = VCF_TO_SNPALN.out.fasta // val(ref+group_meta), fasta
        .combine(CALL_VARIANTS.out.samples, by:0) // val(ref+group_meta), fasta, [meta]
    align_for_tree = align_with_samp_meta
        .filter { it[2].size() >= 3 }
        .map { it[0..1] } // val(ref+group_meta), fasta
    too_few_samp_warnings = align_with_samp_meta // val(ref+group_meta), fasta, [meta]
        .filter { it[2].size() < 3 } 
        .map { [null, it[0].group, it[0].ref, "VARIANT_ANALYSIS", "WARNING", "Not enough samples to build a SNP tree."] } // meta, group_meta, ref_meta, workflow, level, message
    messages = messages.mix(too_few_samp_warnings)
                                                                                
    IQTREE2_SNP ( align_for_tree, [] ) 
    ch_versions = ch_versions.mix(IQTREE2_SNP.out.versions.map{it[0]})

    phylogeny = IQTREE2_SNP.out.phylogeny.map{ [it[0].group, it[0].ref, it[1]] } // group_meta, ref_meta, tree
    snp_align = VCF_TO_SNPALN.out.fasta.map{ [it[0].group, it[0].ref, it[1]] }
    vcf = CALL_VARIANTS.out.vcf.map{ [it[0].group, it[0].ref, it[1]] }
    
    results = CALL_VARIANTS.out.vcf // ref+group_meta, vcf
        .combine(VCF_TO_SNPALN.out.fasta, by:0) // ref+group_meta, vcf, align
        .join(IQTREE2_SNP.out.phylogeny, remainder:true, by:0) // ref+group_meta, vcf, align, tree
        .map{ [it[0].group, it[0].ref] + it[1..3] } // group_meta, ref_meta, vcf, align, tree
        
    emit:
    picard_dict  = REFERENCE_INDEX.out.picard_dict
    samtools_fai = REFERENCE_INDEX.out.samtools_fai
    samtools_gzi = REFERENCE_INDEX.out.samtools_gzi
    bwa_index    = REFERENCE_INDEX.out.bwa_index 
    phylogeny    = phylogeny   // group_meta, ref_meta, tree
    snp_align    = snp_align   // group_meta, ref_meta, fasta
    vcf          = vcf         // group_meta, ref_meta, fasta
    results      = results     // group_meta, ref_meta, vcf, align, tree
    versions     = ch_versions // versions.yml
    messages     = messages    // meta, group_meta, ref_meta, workflow, level, message

}

