process READ2TREE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::read2tree=0.1.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/read2tree:0.1.5--pyhdfd78af_0':
        'biocontainers/read2tree:0.1.5--pyhdfd78af_0' }"

    input:
    tuple val(meta),  path(paired_1, stageAs: 'paired_1_??.fa'), path(paired_2, stageAs: 'paired_2_??.fa'), path(single), path(long_reads), path(markers), path(dna_ref)

    output:
    tuple val(meta), path("${prefix}_read2tree"), emit: out
    tuple val(meta), path("${prefix}_read2tree/tree_merge.nwk"), emit: tree
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def sample_count = paired_1.size() + single.size() + long_reads.size()
    if (sample_count == 1) { // If the process detects only 1 species
        """
        read2tree \\
        ${args} \\
        --threads $task.cpus \\
        --tree \\
        --reads ${paired_1} ${paired_2} ${single} ${long_reads} \\
        --standlone_path ${markers} \\
        --output_path ${prefix}_read2tree

        cat <<-END_VERSIONS > versions.yml
    	"${task.process}":
        	read2tree: \$(echo \$(read2tree --version))
    	END_VERSIONS
        """
    } else { // Otherwise, use multiple species mode
    	"""
    	# This creates the reference folder
    	read2tree --standalone_path ${markers}/ --dna_reference ${dna_ref} --output_path ${prefix}_read2tree --reference

    	# Add each paired end shortread sample
    	for R1 in ${paired_1}; do
        	R2=\$(echo \$R1 | sed 's/^paired_1_/paired_2_/')
        	read2tree \\
            ${args} \\
            --threads $task.cpus \\
        	--standalone_path ${markers}/ \\
            --dna_reference ${dna_ref} \\
        	--output_path ${prefix}_read2tree \\
        	--reads \$R1 \$R2
    	done

    	# Add each single end shortread sample
    	for R1 in ${single}; do
        	read2tree \\
            ${args} \\
            --threads $task.cpus \\
        	--standalone_path ${markers}/ \\
            --dna_reference ${dna_ref} \\
        	--output_path ${prefix}_read2tree \\
        	--reads \$R1
    	done

    	# Add each long read sample
    	for R1 in ${long_reads}; do
        	read2tree \\
            ${args} \\
            --threads $task.cpus \\
        	--standalone_path ${markers}/ \\
            --dna_reference ${dna_ref} \\
        	--output_path ${prefix}_read2tree \\
            --read_type long
        	--reads \$R1
    	done

    	# Build tree
    	read2tree \\
        ${args} \\
        --threads $task.cpus \\
    	--standalone_path ${markers}/ \\
        --dna_reference ${dna_ref} \\
    	--output_path ${prefix}_read2tree -\\
    	-merge_all_mappings \\
    	--tree

    	cat <<-END_VERSIONS > versions.yml
    	"${task.process}":
        	read2tree: \$(echo \$(read2tree --version))
    	END_VERSIONS
    	"""
    }
}
