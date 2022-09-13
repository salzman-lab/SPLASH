process EXTRACT_CBC_UMI {

    tag "${fastq.simpleName}"
    label 'process_very_high'
    conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)

    input:
    tuple val(id), path(fastq)
    val num_fastq_reads

    output:
    tuple val(id), path(outfile)  , emit: seqs

    script:
    outfile                       = "extracted_cbc_umi_${fastq.simpleName}.txt.gz"
    """
    extract_cbc_umi.py \\
        --infile ${fastq} \\
        --id ${id} \\
        --num_fastq_reads ${num_fastq_reads} \\
        --outfile ${outfile}
    """
}
