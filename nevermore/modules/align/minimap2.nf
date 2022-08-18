process minimap2_align {
	label 'align'

	input:
	tuple val(sample), path(reads)
	path(reference)

	output:
	tuple val(sample), path("${sample.id}.sam"), emit: sam

	script:
	def reads = (sample.is_paired) ? "${sample.id}_R1.fastq.gz ${sample.id}_R2.fastq.gz" : "${sample.id}_R1.fastq.gz"
	def mm_options = "--sam-hit-only -t ${task.cpus} -x sr --secondary=yes -a"

	"""
	mkdir -p ${sample.id}/
	minimap2 ${mm_options} --split-prefix ${sample.id}_split ${reference} ${reads} > ${sample.id}/${sample.id}.sam
	"""
}