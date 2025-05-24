python scripts/split_sc_data.py \
    data/Combined_10x_CelSeq2_5cl_data.csv \
    data/Labels.csv \
    data_split_output \
    --test_size 0.25 \
    --random_state 123

Rscript scripts/run_scmap_singleR.R \
    data_split_output/ref_counts.csv \
    data_split_output/ref_labels.csv \
    data_split_output/query_counts.csv \
    data_split_output/query_labels.csv results/

python scripts/run_scanvi.py \
    data_split_output/ref_counts.csv \
    data_split_output/ref_labels.csv \
    data_split_output/query_counts.csv \
    data_split_output/query_labels.csv results/

Rscript scripts/benchmark_analysis.R results/