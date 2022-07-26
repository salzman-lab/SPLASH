
process EXTRACT_CBC_UMI {

    tag "${fastq_id}"
    label 'process_low'
    conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)

    input:
    tuple val(id), path(extracted_R1), path(fastq)
    val num_reads_first_pass

    output:
    tuple val(id), path(outfile)  , emit: seqs

    script:
    outfile                       = "extracted_cbc_umi_${id}.txt.gz"
    """
    extract_cbc_umi.py \\
        --infile ${fastq} \\
        --num_lines ${num_reads_first_pass} \\
        --outfile ${outfile}
    """
}
