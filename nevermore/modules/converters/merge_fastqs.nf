process merge_single_fastqs {
    input:
    tuple val(sample), path(fastqs)

    output:
    tuple val(sample), path("merged/${sample.id}_R1.fastq.gz"), emit: fastq

    script:

    def fastq_in = ""
    def prefix = ""
    if (fastqs instanceof Collection && files.size() == 2) {
        prefix = "cat ${fastqs} |"
    } else {
        fastq_in = "${fastqs[0]}"
    }

    """
    set -e -o pipefail
    mkdir -p merged/

    ${prefix} sortbyname.sh in=${fastq_in} out=merged/${sample.id}_R1.fastq.gz
    """
    // cat *.fastq.gz | sortbyname.sh in=stdin.gz out=merged/${sample.id}_R1.fastq.gz

}
