process collate_feature_counts {
	publishDir "${params.output_dir}", mode: params.publish_mode

	input:
	tuple val(sample), path(count_tables)

	output:
	path("collated/*.txt.gz"), emit: collated, optional: true

	script:
	"""
	mkdir -p collated/
	collate_counts . -o collated/collated -c uniq_scaled
	collate_counts . -o collated/collated -c combined_scaled
	"""
}