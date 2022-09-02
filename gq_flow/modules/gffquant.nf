process run_gffquant_sam {
	publishDir "${params.output_dir}", mode: params.publish_mode

	input:
	tuple val(sample), path(alignments)
	path(gq_db)
	val(gq_params)

	output:
	tuple val(sample), path("${sample}/*.txt.gz"), emit: results

	script:
	def gq_output = "-o ${sample}/${sample}"
	def gq_cmd = "gffquant ${gq_output} ${gq_params} gq_db.sqlite3"

	def mk_aln_sam = ""
	if (samfiles instanceof Collection && samfiles.size() >= 2) {
		mk_aln_sam = "cat ${sample}.sam > alignments.sam \n"
		mk_aln_sam += "grep -v '^@' ${sample}.singles.sam >> alignments.sam"
	} else {
		mk_aln_sam = "ln -s ${alignments[0]} alignments.sam"
	}

	"""
	mkdir -p logs/
	echo 'Copying database...'
	cp -v ${gq_db} gq_db.sqlite3
	echo 'Making alignment stream...'
	${mk_aln_sam}
	cat alignments.sam | ${gq_cmd} - > logs/${sample}.o 2> logs/${sample}.e
	rm -fv gq_db.sqlite3 alignments.sam
	"""
}

process run_gffquant_bam {
	publishDir "${params.output_dir}", mode: params.publish_mode

	input:
	tuple val(sample), path(alignments)
	path(gq_db)
	val(gq_params)

	output:
	tuple val(sample), path("${sample}/*.txt.gz"), emit: results

	script:
	def gq_output = "-o ${sample}/${sample}"
	def gq_cmd = "gffquant --format bam ${gq_output} ${gq_params} gq_db.sqlite3"

	if (params.do_name_sort) {
		gq_cmd = "samtools collate -O ${alignments} -@ ${task.cpus} | ${gq_cmd} -"
	} else {
		gq_cmd = "${gq_cmd} ${alignments}"
	}

	"""
	mkdir -p logs/
	echo 'Copying database...'
	cp -v ${gq_db} gq_db.sqlite3
	${gq_cmd} > logs/${sample}.o 2> logs/${sample}.e
	rm -v gq_db.sqlite3
	"""
}
