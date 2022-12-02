
process ELEMENT_ALIGNMENT {

    tag "${fasta_name}, ${index_name}"
    label 'process_medium'
    conda (params.enable_conda ? 'bioconda::bowtie2=2.4.4 bioconda::samtools=1.15.1 conda-forge::pigz=2.6' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:1744f68fe955578c63054b55309e05b41c37a80d-0' :
        'quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:1744f68fe955578c63054b55309e05b41c37a80d-0' }"

    input:
    tuple path(fasta), val(index)

    output:
    path "*hits*tsv"    , emit: hits

    script:
    fasta_name          = fasta.baseName
    index_name          = file(index).baseName
    outfile             = "${fasta_name}_hits_${index_name}.tsv"
    """
    rm -rf ${outfile}
    echo -e "${fasta_name}\t${fasta_name}_flag_${index_name}\t${fasta_name}_hits_${index_name}\t${fasta_name}_hits_pos_${index_name}" >> ${outfile}

    bowtie2 -f -x ${index} -U ${fasta} -p ${task.cpus} --quiet \\
        | sed '/^@/d' \\
        | cut -f1,2,3,4 \\
        | sort \\
        | uniq \\
        >> ${outfile}
    """
}
