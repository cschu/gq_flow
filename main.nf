#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { fastq_input; bam_input } from "./nevermore/workflows/input"
include { minimap2_align; bwa_mem_align } from "./nevermore/modules/align/minimap2"
include { collate_feature_counts } from "./gq_flow/modules/collate"
include { run_gffquant_sam; run_gffquant_bam } from "./gq_flow/modules/gffquant"

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

if (!params.bam_input_pattern) {
	params.bam_input_pattern = "**.bam"
}
def bam_input_pattern = input_dir + "/" + params.bam_input_pattern

def gq_params = "-m ${params.gq_mode} --ambig_mode ${params.gq_ambig_mode}"
gq_params += (params.gq_strand_specific) ? " --strand_specific" : ""
gq_params += (params.gq_unmarked_orphans) ? " --unmarked_orphans" : ""
gq_params += (params.gq_calc_coverage) ? " --calc_coverage" : ""
gq_params += (params.gq_min_seqlen) ? (" --min_seqlen " + params.gq_min_seqlen) : ""
gq_params += (params.gq_min_identity) ? (" --min_identity " + params.gq_min_identity) : ""

if (params.minimap2_index && params.bwa_mem_index) {
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




process merge_samfiles {
	publishDir params.output_dir, mode: params.publish_mode

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


workflow process_bam_data {
	take:
		aligned_ch
	main:
		run_gffquant_bam(
			aligned_ch,
			params.gq_db,
			gq_params
		)

	emit:
		feature_counts = run_gffquant_bam.out.results
}

workflow process_fastq_data {
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

		run_gffquant_sam(
			aligned_ch,
			params.gq_db,
			gq_params
		)

	emit:
		feature_counts = run_gffquant_sam.out.results

}


workflow {

	feature_count_ch = Channel.empty()

	fastq_input(
		Channel.fromPath(fastq_input_pattern)
	)

	fastq_ch = fastq_input.out.fastqs

	bam_input(
		Channel.fromPath(bam_input_pattern)
	)

	bam_ch = bam_input.out.bamfiles

	process_fastq_data(fastq_ch)

	process_bam_data(bam_ch)

	feature_count_ch = feature_count_ch
		.concat(process_fastq_data.out.feature_counts)
		.concat(process_bam_data.out.feature_counts)
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

	/*
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

	merge_samfiles(aligned_ch)

	align_ch = merge_samfiles.out.alignments

	run_gffquant(
		align_ch,
		params.gq_db,
		gq_params
	)

	feature_count_ch = feature_count_ch
		.concat(run_gffquant.out.results)

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
	*/

}
