## Creating anchor heat maps for top N significant anchors
Use this case if:
1. You have already completed a NOMAD run with annotations(`--run_annotations true`)
2. You want to generate heatmaps for the top 40 anchors, ranked by their anchor significance.
```
nextflow run kaitlinchaung/nomad \
    -r main \
    --latest \
    --input <<samplesheet>> \
    --element_annotations_samplesheet <<element_annotations_samplesheet>> \
    --run_annotations true \
    --results_dir <<Full path to results directory from NOMAD run>> \
    --num_heatmap_anchors 40
```

## Creating anchor heat maps for input list of anchors
Use this case if:
1. You have already completed a NOMAD run with annotations(`--run_annotations true`)
2. You want to generate heatmaps for a list of anchors
    1. This list must be a one column file, with a header titled "anchor"

Please replace params with real paths to files.
```
nextflow run kaitlinchaung/nomad \
    -r main \
    --latest \
    --input <<samplesheet>> \
    --element_annotations_samplesheet <<element_annotations_samplesheet>> \
    --run_annotations true \
    --abundant_stratified_anchors /home/data/results/abundant_stratified_anchors \
    --consensus_fractions /home/data/results/consensus_anchors/*fractions.tab> \
    --anchors_pvals /home/data/results/anchors_pvals.tsv \
    --genome_annotations_anchors /home/data/results/genome_annotations/genome_annotations_anchors.tsv \
    --additional_summary /home/data/results/additional_summary.tsv
```
