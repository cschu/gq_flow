include { run_gffquant } from "../modules/gffquant"
include { minimap2_align; bwa_mem_align } from "../../nevermore/modules/align/minimap2"


workflow process_bam_data {
	take:
		aligned_ch
	main:

		aligned_ch = aligned_ch
			.map { sample, bam ->
				return tuple(sample.id, bam[0])
			}
			.groupTuple(sort: true)

		run_gffquant(
			aligned_ch,
			params.gq_db
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

		run_gffquant(
			aligned_ch,
			params.gq_db
		)

	emit:
		feature_counts = run_gffquant_sam.out.results

}
