process run_gffquant {
	label "gffquant"
	// publishDir "${params.output_dir}", mode: params.publish_mode

	input:
	tuple val(sample), path(alignments)
	path(gq_db)

	output:
	tuple val(sample), path("${sample}/*.txt.gz"), emit: results

	script:
	def gq_output = "-o ${sample}/${sample}"

	def gq_params = "-m ${params.gq_mode} --ambig_mode ${params.gq_ambig_mode}"
	gq_params += (params.gq_strand_specific) ? " --strand_specific" : ""
	gq_params += (params.gq_unmarked_orphans) ? " --unmarked_orphans" : ""
	gq_params += (params.gq_calc_coverage) ? " --calc_coverage" : ""
	gq_params += (params.gq_min_seqlen) ? (" --min_seqlen " + params.gq_min_seqlen) : ""
	gq_params += (params.gq_min_identity) ? (" --min_identity " + params.gq_min_identity) : ""
	gq_params += (params.bam_input_pattern) ? (" --format bam") : " --format sam"

	def gq_cmd = "gffquant ${gq_output} ${gq_params} gq_db.sqlite3"
	
	def mk_aln_sam = ""
	if (params.bam_input_pattern) {

		if (params.do_name_sort) {
			gq_cmd = "samtools collate -@ ${task.cpus} -O ${alignments} tmp/collated_bam | ${gq_cmd} -"
		} else {
			gq_cmd = "${gq_cmd} ${alignments}"
		}

	} else {

		mk_aln_sam += "echo 'Making alignment stream...'\n"
		if (alignments instanceof Collection && alignments.size() >= 2) {
			mk_aln_sam += "cat ${sample}.sam > tmp/alignments.sam \n"
			mk_aln_sam += "grep -v '^@' ${sample}.singles.sam >> tmp/alignments.sam"
		} else {
			mk_aln_sam += "ln -s ../${alignments[0]} tmp/alignments.sam"
		}
		gq_cmd = "cat tmp/alignments.sam | ${gq_cmd} -"

	}

	"""
	set -e -o pipefail
	mkdir -p logs/ tmp/
	echo 'Copying database...'
	cp -v ${gq_db} gq_db.sqlite3
	${mk_aln_sam}
	${gq_cmd} > logs/${sample}.o 2> logs/${sample}.e
	rm -rfv gq_db.sqlite3 tmp/
	"""	

}


process run_gffquant_sam {
	label "gffquant"
	// publishDir "${params.output_dir}", mode: params.publish_mode

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
	if (alignments instanceof Collection && alignments.size() >= 2) {
		mk_aln_sam = "cat ${sample}.sam > alignments.sam \n"
		mk_aln_sam += "grep -v '^@' ${sample}.singles.sam >> alignments.sam"
	} else {
		mk_aln_sam = "ln -s ${alignments[0]} alignments.sam"
	}

	"""
	set -e -o pipefail
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
	label "gffquant"
	// publishDir "${params.output_dir}", mode: params.publish_mode

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
		gq_cmd = "samtools collate -@ ${task.cpus} -O ${alignments} tmp/collated_bam | ${gq_cmd} -"
	} else {
		gq_cmd = "${gq_cmd} ${alignments}"
	}

	"""
	set -e -o pipefail
	mkdir -p logs/ tmp/
	echo 'Copying database...'
	cp -v ${gq_db} gq_db.sqlite3
	${gq_cmd} > logs/${sample}.o 2> logs/${sample}.e
	rm -rf tmp/
	rm -v gq_db.sqlite3
	"""
}
