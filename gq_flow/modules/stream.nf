process stream_minimap2_gffquant {
	// publishDir "${params.output_dir}", mode: params.publish_mode

	input:
	tuple val(sample), path(fastq)
	path(index)
	path(db)
	val(gq_params)

	output:
	tuple val(sample), path("${sample}/*.txt.gz"), emit: results

	script:
	def gq_output = "-o ${sample.id}/${sample.id}"
	def reads = (sample.is_paired) ? "${sample.id}_R1.fastq.gz ${sample.id}_R2.fastq.gz" : "${sample.id}_R1.fastq.gz"
	def mm_options = "--sam-hit-only -t ${task.cpus} -x sr --secondary=yes -a"
	"""
	mkdir -p logs/
	echo 'Copying database...'
	cp -v ${db} gq_db.sqlite3
	minimap2 ${mm_options} --split-prefix ${sample.id}_split ${index} ${reads} | gffquant ${gq_output} ${gq_params} gq_db.sqlite3 - > logs/${sample}.o 2> logs/${sample}.e
	rm -v gq_db.sqlite3
	"""
	// minimap2 --sam-hit-only -t <threads> -x sr --secondary=yes -a [-o <out.sam>] --split-prefix <prefix> <mmi> <reads>
}