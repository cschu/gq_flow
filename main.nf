#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { fastq_input; bam_input } from "./nevermore/workflows/input"
include { minimap2_align } from "./nevermore/modules/align/minimap2"
include { collate_feature_counts } from "./gq_flow/modules/collate"
include { gffquant } from "./gq_flow/modules/gffquant"

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
def bam_input_pattern = input_dir + "/" + "**.bam"

def gq_params = "-m ${params.gq_mode} --ambig_mode ${params.gq_ambig_mode}"
gq_params += (params.gq_strand_specific) ? " --strand_specific" : ""
gq_params += (params.gq_unmarked_orphans) ? " --unmarked_orphans" : ""


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


workflow {

	fastq_input(
		Channel.fromPath(fastq_input_pattern)
	)

	fastq_ch = fastq_input.out.fastqs
	
	feature_count_ch = Channel.empty()

	/*bam_input(
		Channel.fromPath(bam_input_pattern)
	)*/

	minimap2_align(
		fastq_ch,
		params.minimap2_index
	)
	
	aligned_ch = minimap2_align.out.sam
		.map { sample, sam ->
			sample_id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, "")
			return tuple(sample_id, bam)
		}
		.groupTuple(sort: true)

	merge_samfiles(aligned_ch)

	align_ch = merge_samfiles.out.alignments

	run_gffquant(
		align_ch,
		params.minimap2_index,
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

}