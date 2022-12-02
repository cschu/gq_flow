process merge_samfiles {
	// publishDir params.output_dir, mode: params.publish_mode

	input:
		tuple val(sample), path(samfiles)
	output:
		tuple val(sample), path("sam/${sample}.sam"), emit: alignments
	script:
		if (samfiles instanceof Collection && samfiles.size() >= 2) {
			"""
			mkdir -p sam/
			cat ${sample}.sam > sam/${sample}.sam
			grep -v '^@' ${sample}.singles.sam >> sam/${sample}.sam
			"""
		} else {
			"""
			mkdir -p sam/
			ln -s ../${samfiles[0]} sam/${sample}.sam
			"""
		}
}