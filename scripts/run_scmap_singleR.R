#!/usr/bin/env Rscript

# run_scmap_singleR.R
# Usage:
# Rscript run_scmap_singleR.R <ref_counts_file> <ref_labels_file> <query_counts_file> <query_labels_file> <output_dir>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript run_scmap_singleR.R <ref_counts_file> <ref_labels_file> <query_counts_file> <query_labels_file> <output_dir>")
}

ref_counts_file <- args[1]
ref_labels_file <- args[2]
query_counts_file <- args[3]
query_labels_file <- args[4]
output_dir <- args[5]

# Load necessary libraries
library(SingleCellExperiment)
library(scran)
library(Matrix)
library(dplyr)
library(SingleR)
library(BiocParallel)
library(scmap)
library(SummarizedExperiment)

# Initialize a dataframe to store statistics and timings
benchmark_stats_r <- data.frame(
  Dataset = character(),
  Num_Cells = integer(),
  Num_Genes_After_Processing = integer(), # Number of genes might change at different stages
  Num_Cell_Types = integer(),
  Method = character(),
  Execution_Time_Seconds = numeric(),
  stringsAsFactors = FALSE
)

# Set working directory
dir.create(output_dir, showWarnings = FALSE)

# --- 1. Data Preparation ---
message("Loading reference data...")
ref_labels_df <- read.csv(ref_labels_file, header = TRUE, stringsAsFactors = FALSE)
ref_cell_labels_colname <- colnames(ref_labels_df)[1]
ref_cell_labels <- ref_labels_df[[ref_cell_labels_colname]]
ref_counts_df <- read.csv(ref_counts_file, header = TRUE, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)

if (nrow(ref_counts_df) != length(ref_cell_labels)) {
  stop("Mismatch between number of cells in reference counts data and labels.")
}
# Record initial reference set size (number of cells, genes recorded after processing)
benchmark_stats_r <- rbind(benchmark_stats_r, data.frame(
  Dataset = "Reference (Initial)", Num_Cells = nrow(ref_counts_df), Num_Genes_After_Processing = NA, Num_Cell_Types = length(unique(ref_cell_labels)),
  Method = NA, Execution_Time_Seconds = NA
))

ref_counts_matrix <- as(as.matrix(ref_counts_df), "sparseMatrix")
ref_counts_matrix_t <- t(ref_counts_matrix)

message("Creating reference SingleCellExperiment object...")
sce_ref <- SingleCellExperiment(assays = list(counts = ref_counts_matrix_t),
                                colData = DataFrame(cell_id = colnames(ref_counts_matrix_t), label = ref_cell_labels))
rowData(sce_ref)$feature_symbol <- rownames(sce_ref)
colData(sce_ref)$cell_type1 <- as.character(sce_ref$label)

message("Loading query data...")
query_labels_df <- read.csv(query_labels_file, header = TRUE, stringsAsFactors = FALSE)
query_cell_labels_colname <- colnames(query_labels_df)[1]
query_cell_labels <- query_labels_df[[query_cell_labels_colname]]
query_counts_df <- read.csv(query_counts_file, header = TRUE, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)

if (nrow(query_counts_df) != length(query_cell_labels)) {
  stop("Mismatch between number of cells in query counts data and labels.")
}
# Record initial query set size
benchmark_stats_r <- rbind(benchmark_stats_r, data.frame(
  Dataset = "Query (Initial)", Num_Cells = nrow(query_counts_df), Num_Genes_After_Processing = NA, Num_Cell_Types = length(unique(query_cell_labels)),
  Method = NA, Execution_Time_Seconds = NA
))

query_counts_matrix <- as(as.matrix(query_counts_df), "sparseMatrix")
query_counts_matrix_t <- t(query_counts_matrix)

message("Creating query SingleCellExperiment object...")
sce_query <- SingleCellExperiment(assays = list(counts = query_counts_matrix_t),
                                  colData = DataFrame(cell_id = colnames(query_counts_matrix_t), label = query_cell_labels))
rowData(sce_query)$feature_symbol <- rownames(sce_query)

# --- 2. Basic QC and Normalization ---
process_sce <- function(sce, dataset_name) {
  # ... (process_sce function remains unchanged) ...
  message(paste("Processing", dataset_name, "SCE..."))
  message(paste("Initial genes:", nrow(sce), "Initial cells:", ncol(sce)))
  
  message(paste("Filtering low-expressed genes for", dataset_name, "..."))
  ave_counts <- rowMeans(counts(sce))
  keep_genes <- ave_counts > 0.01 
  sce <- sce[keep_genes, ]
  message(paste("Genes after filtering for", dataset_name, ":", nrow(sce)))
  
  message(paste("Normalizing", dataset_name, "data using CPM and log1p..."))
  counts_cpm <- calculateCPM(counts(sce)) 
  logcounts(sce) <- log1p(counts_cpm) 
  
  if (!"logcounts" %in% assayNames(sce)) {
    stop(paste("logcounts assay was not created successfully for", dataset_name))
  }
  message(paste("Normalization complete for", dataset_name, ". logcounts assay added."))
  return(sce)
}

time_start_processing <- Sys.time()
sce_ref_processed <- process_sce(sce_ref, "Reference")
sce_query_processed <- process_sce(sce_query, "Query")
time_end_processing <- Sys.time()
processing_time <- as.numeric(difftime(time_end_processing, time_start_processing, units = "secs"))
benchmark_stats_r <- rbind(benchmark_stats_r, data.frame(
  Dataset = "Combined", Num_Cells = NA, Num_Genes_After_Processing = NA, Num_Cell_Types = NA,
  Method = "Data Processing (Filter & Normalize)", Execution_Time_Seconds = processing_time
))

# Update processed dataset sizes
benchmark_stats_r[benchmark_stats_r$Dataset == "Reference (Initial)", "Num_Genes_After_Processing"] <- nrow(sce_ref_processed)
benchmark_stats_r[benchmark_stats_r$Dataset == "Query (Initial)", "Num_Genes_After_Processing"] <- nrow(sce_query_processed)


# --- 3. Gene Alignment (for SingleR) ---
message("Aligning genes between reference and query sets for SingleR...")
common_genes_for_singleR <- intersect(rownames(sce_ref_processed), rownames(sce_query_processed))
if (length(common_genes_for_singleR) == 0) {
  stop("No common genes found for SingleR.")
}
sce_ref_for_singleR <- sce_ref_processed[common_genes_for_singleR, ]
sce_query_for_singleR <- sce_query_processed[common_genes_for_singleR, ]

message("Reference set for SingleR: ", ncol(sce_ref_for_singleR), " cells, ", nrow(sce_ref_for_singleR), " genes.")
message("Query set for SingleR: ", ncol(sce_query_for_singleR), " cells, ", nrow(sce_query_for_singleR), " genes.")
# Record dataset sizes used by SingleR
benchmark_stats_r <- rbind(benchmark_stats_r, data.frame(
  Dataset = "Reference (for SingleR)", Num_Cells = ncol(sce_ref_for_singleR), Num_Genes_After_Processing = nrow(sce_ref_for_singleR), Num_Cell_Types = length(unique(sce_ref_for_singleR$label)),
  Method = NA, Execution_Time_Seconds = NA
))
benchmark_stats_r <- rbind(benchmark_stats_r, data.frame(
  Dataset = "Query (for SingleR)", Num_Cells = ncol(sce_query_for_singleR), Num_Genes_After_Processing = nrow(sce_query_for_singleR), Num_Cell_Types = length(unique(sce_query_for_singleR$label)),
  Method = NA, Execution_Time_Seconds = NA
))

message("Saving processed query SCE object (after common processing)...")
saveRDS(sce_query_processed, file.path(output_dir, "sce_query_processed.rds"))

# --- 4. SingleR Benchmark ---
message("--- Running SingleR Benchmark ---")
singleR_ref_labels <- sce_ref_for_singleR$label
singleR_ref_logcounts <- assay(sce_ref_for_singleR, "logcounts")
singleR_query_logcounts <- assay(sce_query_for_singleR, "logcounts")
singleR_query_true_labels <- sce_query_for_singleR$label

message("Running SingleR annotation...")
BPPARAM_custom <- BiocParallel::SerialParam()
time_start_singleR <- Sys.time()
singleR_results_obj <- SingleR::SingleR(
  test = singleR_query_logcounts,
  ref = singleR_ref_logcounts,
  labels = singleR_ref_labels,
  genes = "de",
  BPPARAM = BPPARAM_custom
)
time_end_singleR <- Sys.time()
singleR_execution_time <- as.numeric(difftime(time_end_singleR, time_start_singleR, units = "secs"))
benchmark_stats_r <- rbind(benchmark_stats_r, data.frame(
  Dataset = NA, Num_Cells = NA, Num_Genes_After_Processing = NA, Num_Cell_Types = NA,
  Method = "SingleR Annotation", Execution_Time_Seconds = singleR_execution_time
))

message("Extracting SingleR predictions...")
predicted_labels_singleR <- singleR_results_obj$pruned.labels
if (any(is.na(predicted_labels_singleR))) {
  predicted_labels_singleR <- singleR_results_obj$labels
}
if (length(predicted_labels_singleR) != ncol(sce_query_for_singleR)) {
  stop("SingleR prediction length mismatch.")
}
singleR_query_cell_ids <- colnames(sce_query_for_singleR)
results_df_singleR <- data.frame(
  cell_id = singleR_query_cell_ids,
  true_label = as.character(singleR_query_true_labels),
  predicted_label_singleR = as.character(predicted_labels_singleR),
  stringsAsFactors = FALSE
)
write.csv(results_df_singleR, file.path(output_dir, "singleR_predictions.csv"), row.names = FALSE)
message("SingleR predictions saved.")
message("--- SingleR Benchmark Complete ---")

# --- 5. scmap Benchmark (using scmap-cell) ---
message("\n--- Running scmap-cell Benchmark ---")
sce_ref_scmap <- sce_ref_processed
sce_query_scmap <- sce_query_processed

# 5.1 Ensure scmap required columns exist (done or checked during object creation)
colData(sce_ref_scmap)$cell_type1 <- as.character(sce_ref_scmap$label) # Ensure it exists

# 5.2 Feature selection
message("scmap: Selecting features...")
time_start_scmap_feat <- Sys.time()
set.seed(123)
sce_ref_scmap <- selectFeatures(sce_ref_scmap, n_features = 500, suppress_plot = TRUE)
time_end_scmap_feat <- Sys.time()
scmap_feat_time <- as.numeric(difftime(time_end_scmap_feat, time_start_scmap_feat, units = "secs"))
selected_feature_count <- sum(rowData(sce_ref_scmap)$scmap_features, na.rm = TRUE)
message(sprintf("scmap: %d features selected.", selected_feature_count))
if (selected_feature_count == 0) stop("scmap: No features selected.")

# 5.3 Gene alignment and feature transfer
message("scmap: Aligning genes and transferring features...")
common_genes_scmap <- intersect(rownames(sce_ref_scmap), rownames(sce_query_scmap))
if (length(common_genes_scmap) == 0) stop("scmap: No common genes for scmap.")
sce_ref_scmap_common <- sce_ref_scmap[common_genes_scmap, ]
sce_query_scmap_common <- sce_query_scmap[common_genes_scmap, ]
if ("scmap_features" %in% colnames(rowData(sce_ref_scmap_common))) {
  new_row_data_query_scmap <- DataFrame(feature_symbol = rownames(sce_query_scmap_common))
  new_row_data_query_scmap$scmap_features <- rowData(sce_ref_scmap_common)$scmap_features
  rowData(sce_query_scmap_common) <- new_row_data_query_scmap
} else {
  warning("scmap: 'scmap_features' not found after common gene selection.")
}
# Record dataset sizes used by scmap
benchmark_stats_r <- rbind(benchmark_stats_r, data.frame(
  Dataset = "Reference (for scmap)", Num_Cells = ncol(sce_ref_scmap_common), Num_Genes_After_Processing = selected_feature_count, Num_Cell_Types = length(unique(sce_ref_scmap_common$label)),
  Method = NA, Execution_Time_Seconds = NA
))
benchmark_stats_r <- rbind(benchmark_stats_r, data.frame(
  Dataset = "Query (for scmap)", Num_Cells = ncol(sce_query_scmap_common), Num_Genes_After_Processing = nrow(sce_query_scmap_common), Num_Cell_Types = length(unique(sce_query_scmap_common$label)), # Genes here are common, not just selected
  Method = NA, Execution_Time_Seconds = NA
))


# 5.4 Build scmap-cell index
message("scmap: Building scmap-cell index...")
time_start_scmap_index <- Sys.time()
set.seed(1)
sce_ref_scmap_indexed <- indexCell(sce_ref_scmap_common, k = 10)
time_end_scmap_index <- Sys.time()
scmap_index_time <- as.numeric(difftime(time_end_scmap_index, time_start_scmap_index, units = "secs"))
if (!"scmap_cell_index" %in% names(metadata(sce_ref_scmap_indexed))) stop("scmap: Indexing failed.")
scmap_cell_index_data <- metadata(sce_ref_scmap_indexed)$scmap_cell_index

# 5.5 Project and annotate using scmap-cell
message("scmap: Projecting query data...")
time_start_scmap_predict <- Sys.time()
scmapCell_results_obj <- scmapCell(
  projection = sce_query_scmap_common,
  index_list = list(scmap_ref = scmap_cell_index_data)
)
reference_cell_types_for_scmap <- as.character(colData(sce_ref_scmap_indexed)$cell_type1)
scmapCell_clusters_obj <- scmapCell2Cluster(
  scmapCell_results_obj,
  list(scmap_ref = reference_cell_types_for_scmap), w = 1
)
time_end_scmap_predict <- Sys.time()
scmap_predict_time <- as.numeric(difftime(time_end_scmap_predict, time_start_scmap_predict, units = "secs"))

# Combine scmap times
scmap_total_time <- scmap_feat_time + scmap_index_time + scmap_predict_time
benchmark_stats_r <- rbind(benchmark_stats_r, data.frame(
  Dataset = NA, Num_Cells = NA, Num_Genes_After_Processing = NA, Num_Cell_Types = NA,
  Method = "scmap-cell (FeatSel+Index+Predict)", Execution_Time_Seconds = scmap_total_time
))
benchmark_stats_r <- rbind(benchmark_stats_r, data.frame(
  Dataset = NA, Num_Cells = NA, Num_Genes_After_Processing = NA, Num_Cell_Types = NA,
  Method = "scmap-cell (Feature Selection)", Execution_Time_Seconds = scmap_feat_time
))
benchmark_stats_r <- rbind(benchmark_stats_r, data.frame(
  Dataset = NA, Num_Cells = NA, Num_Genes_After_Processing = NA, Num_Cell_Types = NA,
  Method = "scmap-cell (Indexing)", Execution_Time_Seconds = scmap_index_time
))
benchmark_stats_r <- rbind(benchmark_stats_r, data.frame(
  Dataset = NA, Num_Cells = NA, Num_Genes_After_Processing = NA, Num_Cell_Types = NA,
  Method = "scmap-cell (Prediction)", Execution_Time_Seconds = scmap_predict_time
))


# 5.7 Extract prediction results
message("scmap: Extracting predictions...")
if (!"combined_labs" %in% names(scmapCell_clusters_obj)) {
  if ("scmap_cluster_labs" %in% names(scmapCell_clusters_obj) && "scmap_ref" %in% names(scmapCell_clusters_obj$scmap_cluster_labs)) {
    predicted_labels_scmap_cell <- as.character(scmapCell_clusters_obj$scmap_cluster_labs$scmap_ref)
  } else { stop("scmap: Could not find predictions.") }
} else {
  predicted_labels_scmap_cell <- as.character(scmapCell_clusters_obj$combined_labs)
}
if (length(predicted_labels_scmap_cell) != ncol(sce_query_scmap_common)) {
  stop("scmap: Prediction length mismatch.")
}

# 5.8 Save prediction results
scmap_query_cell_ids <- colnames(sce_query_scmap_common)
scmap_query_true_labels_aligned <- as.character(colData(sce_query_scmap_common)$label)
results_df_scmap_cell <- data.frame(
  cell_id = scmap_query_cell_ids,
  true_label = scmap_query_true_labels_aligned,
  predicted_label_scmap = predicted_labels_scmap_cell,
  stringsAsFactors = FALSE
)
write.csv(results_df_scmap_cell, file.path(output_dir, "scmap_predictions.csv"), row.names = FALSE)
message("scmap-cell predictions saved.")
message("--- scmap-cell Benchmark Complete ---")

# --- 6. Save Statistics ---
# Prepare final benchmark statistics: include method timings and dataset size descriptors.
benchmark_stats_r_final <- benchmark_stats_r[!is.na(benchmark_stats_r$Method) | 
                                               (benchmark_stats_r$Dataset %in% c("Reference (Initial)", "Query (Initial)", 
                                                                                 "Reference (for SingleR)", "Query (for SingleR)", 
                                                                                 "Reference (for scmap)", "Query (for scmap)")), ]

write.csv(benchmark_stats_r_final, file.path(output_dir, "r_benchmark_stats.csv"), row.names = FALSE)
message("R benchmark statistics and timings saved to ", file.path(output_dir, "r_benchmark_stats.csv"))

message("\nAll R benchmarks complete.")
