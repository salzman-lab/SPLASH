
process ANCHOR_HEATMAPS {

    publishDir "${params.outdir}/anchor_heatmaps", mode: 'copy'
    tag "heatmaps"
    label 'process_high'
    conda (params.enable_conda ? "conda-forge::python=3.9.5 numpy conda-forge::matplotlib pandas scipy seaborn" : null)

    input:
    val use_heatmap_anchor_list
    path heatmap_anchor_list
    path abundant_stratified_anchors
    path consensus_fractions
    val dataset
    val kmer_size
    path anchors_pvals
    path anchor_Cjs
    val num_heatmap_anchors
    path samplesheet
    path additional_summary
    path genome_annotations_anchors
    path results_dir

    output:
    path "*.png"                        , emit: png                       , optional: true
    path "*.tsv"                        , emit: tsv                       , optional: true

    script:
    def heatmap_anchor_list             = (use_heatmap_anchor_list==true) ? "--heatmap_anchor_list ${heatmap_anchor_list}" : ""
    def use_additional_summary          = (additional_summary)            ? "--additional_summary ${additional_summary}" : ""
    def use_genome_annotations_anchors  = (genome_annotations_anchors)    ? "--genome_annotations_anchors ${genome_annotations_anchors}" : ""

    outfile_contingency_table           = "heatmap_anchors_contingency_table.csv"
    outfile_skipped_anchors             = "heatmap_skipped_anchors.tsv"


    """
    anchor_heatmaps.py \\
        --dataset ${dataset} \\
        --kmer_size ${kmer_size} \\
        --anchor_pvals ${anchors_pvals} \\
        --anchor_Cjs ${anchor_Cjs} \\
        --num_heatmap_anchors ${num_heatmap_anchors} \\
        --samplesheet ${samplesheet} \\
        --outfile_contingency_table ${outfile_contingency_table} \\
        --outfile_skipped_anchors ${outfile_skipped_anchors} \\
        ${heatmap_anchor_list} \\
        ${use_additional_summary} \\
        ${use_genome_annotations_anchors} \\
    """
}
