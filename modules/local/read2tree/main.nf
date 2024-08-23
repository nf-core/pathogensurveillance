process READ2TREE {
    tag "$meta.id"
    label 'process_medium'
    errorStrategy 'ignore'

    conda "bioconda::read2tree=0.1.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/read2tree:0.1.5--pyhdfd78af_0':
        'biocontainers/read2tree:0.1.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), val(pair_meta), path(paired_1), path(paired_2), val(single_meta), path(single), val(long_meta), path(long_reads), path(markers), path(dna_ref)

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
        FORWARD=(${paired_1})
        REVERSE=(${paired_2})
        IDS=(${pair_meta.collect{it.id}.join(' ')})
    	for i in \${!FORWARD[@]}; do
        	read2tree \\
            ${args} \\
            --threads $task.cpus \\
        	--standalone_path ${markers}/ \\
            --dna_reference ${dna_ref} \\
        	--output_path ${prefix}_read2tree \\
        	--reads \${FORWARD[\$i]} \${REVERSE[\$i]} \\
            --species_name \${IDS[\$i]}
    	done

    	# Add each single end shortread sample
        SINGLE=(${single})
        IDS=(${single_meta.collect{it.id}.join(' ')})
    	for i in \${!SINGLE[@]}; do
        	read2tree \\
            ${args} \\
            --threads $task.cpus \\
        	--standalone_path ${markers}/ \\
            --dna_reference ${dna_ref} \\
        	--output_path ${prefix}_read2tree \\
        	--reads \${SINGLE[\$i]} \\
            --species_name \${IDS[\$i]}
    	done

    	# Add each long read sample
        LONG=(${long_reads})
        IDS=(${long_meta.collect{it.id}.join(' ')})
    	for i in \${!LONG[@]}; do
        	read2tree \\
            ${args} \\
            --threads $task.cpus \\
        	--standalone_path ${markers}/ \\
            --dna_reference ${dna_ref} \\
        	--output_path ${prefix}_read2tree \\
            --read_type long \\
        	--reads \${LONG[\$i]} \\
            --species_name \${IDS[\$i]}
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
