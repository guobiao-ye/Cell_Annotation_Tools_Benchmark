# Cell_Annotation_Tools_Benchmark

Welcome to the **Cell_Annotation_Tools_Benchmark** repository, a comprehensive benchmark study evaluating the performance of three automated cell type annotation tools—**SingleR**, **scmap**, and **scANVI**—for single-cell RNA-sequencing (scRNA-seq) data. This project systematically assesses these tools under varying reference data availability, batch heterogeneity, and computational constraints, providing practical guidance for tool selection in real-world scRNA-seq analysis.

This repository contains the code, data processing pipeline, and visualization scripts used to generate the results presented in the associated research report. The study focuses on pragmatic robustness, addressing challenges such as scarce reference data, batch effects, and computational scalability.

## Table of Contents
1. [Project Overview](#project-overview)
2. [Repository Structure](#repository-structure)
3. [Dependencies](#dependencies)
4. [Installation](#installation)
5. [Usage](#usage)
6. [Data](#data)
7. [Pipeline Details](#pipeline-details)
8. [Key Findings](#key-findings)
9. [Visualization](#visualization)
10. [Contributing](#contributing)
11. [License](#license)
12. [References](#references)

## Project Overview
Single-cell RNA-sequencing (scRNA-seq) enables high-resolution analysis of cellular heterogeneity, but accurate cell type annotation remains a critical challenge. Automated annotation tools streamline this process, yet their performance varies under real-world constraints such as limited reference data, batch effects, and computational resources. This project benchmarks three representative tools:

- **SingleR (v2.4.1)**: A correlation-based method using Spearman rank correlation on highly variable genes.
- **scmap (v1.24.0, scmap-cell)**: A projection-based method employing k-nearest neighbors or centroids with cosine similarity.
- **scANVI (v1.3.1)**: A deep learning method using a variational autoencoder (VAE) for probabilistic data integration and batch correction.

The benchmark evaluates these tools across:
- **Reference Scale**: From extremely scarce (1.25% of cells) to abundant (75% of cells) reference sets.
- **Batch Heterogeneity**: Using a single CEL-Seq2 dataset and a combined 10x Chromium + CEL-Seq2 dataset to simulate batch effects.
- **Computational Efficiency**: Measuring execution time, CPU RAM, and GPU VRAM (for scANVI).

The results provide a pragmatic decision framework for selecting the most suitable tool based on data characteristics and computational constraints.

## Project Structure
```
Cell_Annotation_Tools_Benchmark/
├── data/                           # Input data (not included, download from GEO)
├── data_split_output/             # Output directory for split datasets
├── results/                       # Output directory for predictions and benchmark results
├── scripts/                       # Core pipeline scripts
│   ├── split_sc_data.py           # Splits dataset into reference and query sets
│   ├── run_scmap_singleR.R        # Runs SingleR and scmap annotations
│   ├── run_scanvi.py              # Runs scANVI annotation
│   ├── benchmark_analysis.R       # Analyzes predictions and generates metrics
├── visualization/                 # Visualization scripts
│   ├── celseq2_accuracy_comparison.R           # Accuracy plot for CEL-Seq2 dataset
│   ├── combined_accuracy_comparison.R          # Accuracy plot for combined dataset
│   ├── combined_time_comparison.R             # Execution time comparison
│   ├── combined_stacked_time_stages.R         # Stacked time breakdown by stage
│   ├── combined_stacked_CPU_GPU.R             # Resource consumption plot
├── pipeline.sh                    # Main pipeline script to run all steps
├── README.md                      # This file
```

## Dependencies
To run the pipeline, ensure the following software and packages are installed:

### Software
- **Python** (v3.12 or later)
- **R** (v4.3.2 or later)

### Python Packages
```bash
pip install pandas numpy scikit-learn scanpy scvi-tools torch
```

### R Packages
```R
install.packages(c("dplyr", "ggplot2", "MLmetrics", "vcd", "patchwork", "cowplot", "tidyr", "scales"))
BiocManager::install(c("SingleCellExperiment", "scran", "Matrix", "SingleR", "BiocParallel", "scmap", "SummarizedExperiment", "umap"))
```

### Hardware
- For **scANVI**, a GPU is recommended to reduce training time, though CPU execution is supported.

## Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/guobiao-ve/Cell_Annotation_Tools_Benchmark.git
   cd Cell_Annotation_Tools_Benchmark
   ```

2. Install dependencies:
   - Install Python packages using the command above.
   - Install R packages using the provided R commands.

3. Download the datasets:
   - **CellBench 10X**: GEO Accession [GSM3618014](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3618014)
   - **CellBench CEL-Seq2**: GEO Accessions [GSM3618022](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3618022), [GSM3618023](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3618023), [GSM3618024](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3618024)
   - Place the processed data in the `data/` directory as `Combined_10x_CelSeq2_5cl_data.csv` (expression matrix) and `Labels.csv` (cell type labels).

## Usage
1. **Prepare Input Data**:
   - Ensure `data/Combined_10x_CelSeq2_5cl_data.csv` (cells x genes, cell IDs as row names) and `data/Labels.csv` (cell labels, column name "x") are in place.
   - The pipeline assumes these files are preprocessed and aligned.

2. **Run the Pipeline**:
   Execute the main pipeline script to split data, run annotations, and analyze results:
   ```bash
   bash pipeline.sh
   ```
   This script:
   - Splits the dataset into reference and query sets (`scripts/split_sc_data.py`).
   - Runs SingleR and scmap annotations (`scripts/run_scmap_singleR.R`).
   - Runs scANVI annotation (`scripts/run_scanvi.py`).
   - Analyzes predictions and generates metrics (`scripts/benchmark_analysis.R`).

3. **Generate Visualizations**:
   Run the visualization scripts in the `visualization/` directory to produce plots:
   ```R
   Rscript visualization/celseq2_accuracy_comparison.R
   Rscript visualization/combined_accuracy_comparison.R
   Rscript visualization/combined_time_comparison.R
   Rscript visualization/combined_stacked_time_stages.R
   Rscript visualization/combined_stacked_CPU_GPU.R
   ```
   Outputs are saved as PNG files in the `visualization/` directory.

4. **Output Files**:
   - **Data Splits**: Stored in `data_split_output/` (e.g., `ref_counts.csv`, `query_labels.csv`).
   - **Predictions**: Stored in `results/` (e.g., `singleR_predictions.csv`, `scanvi_predictions.csv`, `scmap_predictions.csv`).
   - **Metrics**: Stored in `results/` (e.g., `overall_metrics.csv`, `class_metrics.csv`, `benchmark_report.txt`).
   - **Visualizations**: Stored in `visualization/` (e.g., `combined_accuracy_comparison.png`).

## Data
The benchmark uses two public scRNA-seq datasets:
- **CellBench 10X** (GSM3618014): 3,803 cells, 5 cell types (A549, H1975, HCC827, etc.), 10x Chromium platform.
- **CellBench CEL-Seq2** (GSM3618022-24): 570 cells, 5 cell types, CEL-Seq2 platform.
- **Combined Dataset**: Merges 10x and CEL-Seq2 (4,373 cells) to simulate batch effects.

The pipeline expects a combined expression matrix (`Combined_10x_CelSeq2_5cl_data.csv`) and a label file (`Labels.csv`) with consistent cell ID ordering.

## Pipeline Details
The pipeline consists of four main scripts, orchestrated by `pipeline.sh`:

1. **`split_sc_data.py`**:
   - Splits the input dataset into reference and query sets using stratified sampling (`test_size=0.25`, `random_state=123`).
   - Ensures cell ID and label alignment.
   - Outputs: `ref_counts.csv`, `ref_labels.csv`, `query_counts.csv`, `query_labels.csv`.

2. **`run_scmap_singleR.R`**:
   - Processes reference and query data (gene filtering: mean counts > 0.01, log-normalization).
   - Runs SingleR with default parameters and scmap with 500 features and k=10 neighbors.
   - Saves predictions and processed query data (`sce_query_processed.rds`).

3. **`run_scanvi.py`**:
   - Processes data using scanpy (normalization: target sum 1e4, log1p, 2000 highly variable genes).
   - Trains an scVI model (200 epochs, 2 layers, 30 latent dimensions), converts to scANVI, and fine-tunes (100 epochs).
   - Saves predictions and processed query data (`adata_query_processed.h5ad`).

4. **`benchmark_analysis.R`**:
   - Merges predictions for common cells.
   - Computes metrics: accuracy, Cohen's Kappa, per-class precision, recall, F1-score.
   - Generates visualizations: confusion matrices, UMAP plots, performance bar plots.
   - Outputs: `overall_metrics.csv`, `class_metrics.csv`, `benchmark_report.txt`, and PNG files.

## Key Findings
The benchmark provides insights into tool performance under varying conditions:

### Reference Data Scarcity (CEL-Seq2 Dataset)
- **SingleR**: Highly robust at low reference proportions (e.g., 2.5%, ~14 cells, accuracy > 0.98) due to rank-based correlation.
- **scANVI**: Struggles with sparse references (accuracy 0.624–0.811 at 2.5–5%) due to VAE overfitting.
- **scmap**: Poor performance with sparse references (accuracy 0.309–0.601 at 2.5–5%) due to unstable feature selection.

### Batch Heterogeneity (Combined Dataset)
- **scANVI**: Excels with large, heterogeneous references (accuracy = 1 at 75%, ~3280 cells) due to batch correction via VAE.
- **SingleR**: Robust (accuracy 0.956–0.981 across proportions) but requires preprocessing for batch effects.
- **scmap**: High accuracy with large references (0.999 at 75%) but sensitive to batch effects at low proportions.

### Computational Efficiency
- **SingleR**: Fastest (1.5–7.5 seconds, 0.6–3.2 GB CPU RAM), ideal for quick annotations.
- **scmap**: Fast for small queries but scales poorly with large query sets (up to 11 seconds).
- **scANVI**: High computational cost (up to 44.3 seconds, 12.6 GB CPU RAM, 6.7 GB GPU VRAM) but justified for complex datasets.

### Cross-Batch Generalization
- All tools achieve >98% accuracy when trained on 10x and tested on CEL-Seq2, indicating robustness to batch differences when reference data captures cell type signatures.

### Decision Framework
- **Sparse, Homogeneous Reference**: Use **SingleR** for speed and accuracy.
- **Large, Heterogeneous Reference**: Use **scANVI** for superior batch correction.
- **Moderate Reference, Large Queries**: Use **scmap** but monitor prediction scalability.

## Visualization
The `visualization/` directory contains scripts to generate key figures:
- **Accuracy Comparison**: Bar plots for CEL-Seq2 (`celseq2_accuracy_comparison.R`) and combined datasets (`combined_accuracy_comparison.R`).
- **Execution Time**: Line plot comparing total execution time (`combined_time_comparison.R`).
- **Time Breakdown**: Stacked bar plot showing time per stage (`combined_stacked_time_stages.R`).
- **Resource Consumption**: Stacked bar plot for CPU RAM and GPU VRAM (`combined_stacked_CPU_GPU.R`).

Run these scripts to reproduce figures as shown in the research report.

## Contributing
Contributions are welcome! Please:
1. Fork the repository.
2. Create a feature branch (`git checkout -b feature/new-analysis`).
3. Commit changes (`git commit -m 'Add new analysis'`).
4. Push to the branch (`git push origin feature/new-analysis`).
5. Open a pull request.

For issues or suggestions, please use the [Issues](https://github.com/guobiao-ve/Cell_Annotation_Tools_Benchmark/issues) tab.

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## References
- Abdelaal, T., et al. (2019). A comparison of automatic cell identification methods for single-cell RNA sequencing data. *Genome Biology*, 20(1), 194.
- Aran, D., et al. (2019). Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. *Nature Immunology*, 20(2), 163–172.
- Kiselev, V. Y., et al. (2018). scmap: projection of single-cell RNA-seq data across data sets. *Nature Methods*, 15(5), 359–362.
- Tian, L., et al. (2019). Benchmarking single cell RNA-sequencing analysis pipelines using mixture control experiments. *Nature Methods*, 16(6), 479–487.
- Xu, C., et al. (2021). Probabilistic harmonization and annotation of single-cell transcriptomics data with deep generative models. *Molecular Systems Biology*, 17(1), e9620.

---

Thank you for exploring the **Cell_Annotation_Tools_Benchmark** repository! For questions, contact the me via GitHub.
