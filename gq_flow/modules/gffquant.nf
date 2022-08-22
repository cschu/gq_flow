process run_gffquant {
	publishDir "${params.output_dir}", mode: params.publish_mode

	input:
	tuple val(sample), path(bam)
	path(db)
	val(gq_params)

	output:
	tuple val(sample), path("${sample}/*.txt.gz"), emit: results

	script:
	def gq_output = "-o ${sample}/${sample}"
	"""
	mkdir -p logs/
	echo 'Copying database...'
	cp -v ${db} gq_db.sqlite3
	gffquant ${gq_output} ${gq_params} gq_db.sqlite3 - > logs/${sample}.o 2> logs/${sample}.e
	rm -v gq_db.sqlite3
	"""
}
