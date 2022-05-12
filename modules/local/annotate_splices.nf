
process ANNOTATE_SPLICES {

    label 'process_medium'
    conda (params.enable_conda ? 'bioconda::star=2.7.9a' : null)

    input:
    tuple val(fasta_name), path(junctions)
    path exon_pickle
    path splice_pickle

    output:
    path "*tsv" , emit: tsv

    script:
    outfile     = "ann_splices_${fasta_name}.tsv"
    """
    awk -v OFS='\t' '{print \$1, \$2-1, \$3+1, \$4, \$5, \$6, \$7, \$8, \$9}' ${junctions} \
        > positionModified_${fasta_name}.txt

    ann_splices.py \
        --in_file positionModified_${fasta_name}.txt \
        --out_file ${outfile} \
        --exon_pickle ${exon_pickle} \
        --splice_pickle ${splice_pickle}
    """
}
