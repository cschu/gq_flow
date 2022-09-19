#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { fastq_input; bam_input } from "./nevermore/workflows/input"
// include { minimap2_align; bwa_mem_align } from "./nevermore/modules/align/minimap2"
include { collate_feature_counts } from "./gq_flow/modules/collate"
// include { run_gffquant_sam; run_gffquant_bam } from "./gq_flow/modules/gffquant"

if (params.input_dir && params.remote_input_dir) {
	log.info """
		Cannot process both --input_dir and --remote_input_dir. Please check input parameters.
	""".stripIndent()
	exit 1
} else if (!params.input_dir && !params.remote_input_dir) {
	log.info """
		Neither --input_dir nor --remote_input_dir set.
	""".stripIndent()
	exit 1
}

def input_dir = (params.input_dir) ? params.input_dir : params.remote_input_dir
def fastq_input_pattern = input_dir + "/" + "**[._]{fastq.gz,fq.gz}"

def bam_input_pattern = null
if (params.bam_input_pattern) {
	bam_input_pattern = input_dir + "/" + params.bam_input_pattern
} else if (params.minimap2_index && params.bwa_mem_index) {
	log.info """
	Please only specify one index (--bwa_mem_index or --minimap2_index).
	""".stripIndent()
	exit 1
} else if (!params.minimap2_index && !params.bwa_mem_index) {
	log.info """
	Neither --bwa_mem_index nor --minimap2_index specified.
	""".stripIndent()
	exit 1
}


workflow {

	feature_count_ch = Channel.empty()

	if (bam_input_pattern) {

		bam_input(
			Channel.fromPath(bam_input_pattern)
		)

		bam_ch = bam_input.out.bamfiles

		process_bam_data(bam_ch)

		feature_count_ch = feature_count_ch
			.concat(process_bam_data.out.feature_counts)

	} else if (params.test_transfer || params.minimap2_index || params.bwa_mem_index) {
		fastq_input(
			Channel.fromPath(fastq_input_pattern)
		)

		fastq_ch = fastq_input.out.fastqs

		if (!params.test_transfer) {
			process_fastq_data(fastq_ch)
		
		feature_count_ch = feature_count_ch
			.concat(process_fastq_data.out.feature_counts)
		}
	}

	if (!params.test_transfer) {
		feature_count_ch = feature_count_ch
			.map { sample, files -> return files }
			.flatten()
			.map { file ->
				def category = file.name
					.replaceAll(/\.txt$/, "")
					.replaceAll(/.+\./, "")
				return tuple(category, file)
			}
			.groupTuple(sort: true)

		if (!params.no_collate) {
			collate_feature_counts(feature_count_ch)
		}
	}
}
