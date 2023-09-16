#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { fastq_input; bam_input } from "./nevermore/workflows/input"
include { run_gffquant; collate_feature_counts } from "./nevermore/modules/profilers/gffquant"
include { minimap2_align; bwa_mem_align } from "./nevermore/modules/align/sam_align"
include { db_filter; db2bed3; readcount } from "./nevermore/modules/align/helpers"


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


workflow align_reads {
	take:
		fastq_ch
	main:
		aligned_ch = Channel.empty()

		if (params.minimap2_index) {

			minimap2_align(
				fastq_ch,
				params.minimap2_index
			)

			aligned_ch = minimap2_align.out.sam

		} else if (params.bwa_mem_index) {

			bwa_mem_align(
				fastq_ch,
				params.bwa_mem_index
			)

			aligned_ch = bwa_mem_align.out.sam

		}

		aligned_ch = aligned_ch
			.map { sample, sam ->
				sample_id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, "")
				return tuple(sample_id, sam)
			}
			.groupTuple(sort: true)

	emit:
		alignments = aligned_ch	

}


workflow {

	alignment_ch = Channel.empty()
	feature_count_ch = Channel.empty()

	if (params.minimap2_index || params.bwa_mem_index) {

		// fastq_input(
		// 	Channel.fromPath(fastq_input_pattern)
		// )
		fastq_input(
			Channel.fromPath(input_dir + "/*", type: "dir")
		)	
		fastq_ch = fastq_input.out.fastqs

		align_reads(fastq_ch)
		alignment_ch = align_reads.out.alignments

	} else if (bam_input_pattern) {

		bam_input(
			Channel.fromPath(bam_input_pattern)
		)
		alignment_ch = bam_input.out.bamfiles
			.map { sample, bam ->
				return tuple(sample.id, bam[0])
			}
			.groupTuple(sort: true)
	}

	// readcount(alignment_ch)

	// alignment_ch = alignment_ch
	// 	.join(readcount.out.readcounts)

	run_gffquant(
		alignment_ch,
		params.gq_db
	)

	feature_count_ch = run_gffquant.out.results
		.map { sample, files -> return files }
		.flatten()
		.filter { !it.name.endsWith("Counter.txt.gz") }
		.filter { params.collate_gene_counts || !it.name.endsWith("gene_counts.txt.gz") }
		.map { file ->
			def category = file.name
				.replaceAll(/\.txt\.gz$/, "")
				.replaceAll(/.+\./, "")
			return tuple(category, file)
		}
		.groupTuple(sort: true)
		.combine(
			Channel.from(params.gq_collate_columns.split(","))
		)

	feature_count_ch.view()

	if (!params.no_collate) {
		collate_feature_counts(feature_count_ch)
	}
	
}
