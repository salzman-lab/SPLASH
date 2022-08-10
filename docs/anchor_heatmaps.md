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

The list of anchors file (`anchor_list.tsv`) might look like this:
```
anchor
AAAAAAAAAAAA
CCCCCCCCCCCC
```

Please replace params with real paths to files.
```
nextflow run kaitlinchaung/nomad \
    -r main \
    --latest \
    --input <<samplesheet>> \
    --element_annotations_samplesheet <<element_annotations_samplesheet>> \
    --run_annotations true \
    --results_dir <<Full path to results directory from NOMAD run>> \
    --use_heatmap_anchor_list true \
    --run_anchor_heatmaps true \
    --heatmap_anchor_list anchor_list.tsv

```

## Outputs
Outputs would be found in `./results/anchor_heatmaps`.
