process EXTRACT_CBC_UMI {

    tag "${id}"
    label 'process_low'
    conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)

    input:
    tuple val(id), path(fastq)
    val num_fastq_reads

    output:
    tuple val(id), path(outfile)  , emit: seqs

    script:
    outfile                       = "extracted_cbc_umi_${id}.txt.gz"
    """
    extract_cbc_umi.py \\
        --infile ${fastq} \\
        --num_lines ${num_fastq_reads} \\
        --id ${id} \\
        --outfile ${outfile}
    """
}
