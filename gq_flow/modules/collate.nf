params.gq_collate_columns = "uniq_scaled,combined_scaled"

process collate_feature_counts {
	// publishDir "${params.output_dir}", mode: params.publish_mode

	input:
	tuple val(sample), path(count_tables)

	output:
	path("collated/*.txt.gz"), emit: collated, optional: true

	script:
	"""
	mkdir -p collated/

	for column in \$(echo ${params.gq_collate_columns} | sed 's/,/ /g'); do
		collate_counts . -o collated/collated -c \$column
	done
	"""
	// collate_counts . -o collated/collated -c uniq_scaled
	// collate_counts . -o collated/collated -c combined_scaled
}