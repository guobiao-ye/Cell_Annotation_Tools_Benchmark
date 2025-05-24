#!/usr/bin/env Rscript

# benchmark_analysis.R
# Usage:
# Rscript benchmark_analysis.R <results_dir>
# Example:
# Rscript benchmark_analysis.R results/  (assuming output from run_scmap_singleR.R is in results/)
# or specify the directory containing all prediction CSV files

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript benchmark_analysis.R <results_dir_containing_predictions_and_sce_query>")
}

results_dir <- args[1]

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(MLmetrics)
library(vcd) # For Kappa
library(SingleCellExperiment) # For UMAP on sce_query_processed.rds
library(umap) # umap package
library(patchwork)
library(cowplot) # For plot_grid, if needed (patchwork is often preferred)
library(tidyr)   # For pivot_wider

# --- 1. Load prediction results ---
message("Loading prediction results...")
singleR_file <- file.path(results_dir, "singleR_predictions.csv")
scanvi_file <- file.path(results_dir, "scanvi_predictions.csv") # assuming scANVI results are also in the same directory or accessible
scmap_file <- file.path(results_dir, "scmap_predictions.csv")

if (!file.exists(singleR_file)) stop(paste("SingleR predictions file not found:", singleR_file))
if (!file.exists(scanvi_file)) stop(paste("scANVI predictions file not found:", scanvi_file))
if (!file.exists(scmap_file)) stop(paste("scmap predictions file not found:", scmap_file))

singleR_results_raw <- read.csv(singleR_file, stringsAsFactors = FALSE)
scanvi_results_raw <- read.csv(scanvi_file, stringsAsFactors = FALSE)
scmap_results_raw <- read.csv(scmap_file, stringsAsFactors = FALSE)

# Rename columns for consistency (assuming first column is cell_id, second is true_label)
# and the prediction column is named predicted_label_METHOD
standardize_colnames <- function(df, method_name_pred_col) {
  # Ensure cell_id, true_label, and prediction columns exist and are named correctly
  if (!"cell_id" %in% colnames(df)) {
    if (ncol(df) > 0) {
      message(paste("Assuming first column of", deparse(substitute(df)), "is 'cell_id'. Renaming."))
      colnames(df)[1] <- "cell_id"
    } else {
      stop(paste("DataFrame", deparse(substitute(df)), "is empty or has no columns."))
    }
  }
  if (!"true_label" %in% colnames(df)) {
    if (ncol(df) > 1) {
      message(paste("Assuming second column of", deparse(substitute(df)), "is 'true_label'. Renaming."))
      colnames(df)[2] <- "true_label"
    } else {
      stop(paste("DataFrame", deparse(substitute(df)), "does not have enough columns for 'true_label'."))
    }
  }
  # Find the column containing "predicted_label_" and ensure it's the correct prediction column
  pred_col_idx <- grep(paste0("predicted_label_", method_name_pred_col), colnames(df), ignore.case = TRUE)
  if (length(pred_col_idx) == 1) {
    colnames(df)[pred_col_idx] <- paste0("predicted_label_", method_name_pred_col)
  } else if (method_name_pred_col == "scmap" && "predicted_label_scmap" %in% colnames(df)) {
    # For scmap_predictions.csv which might already have the correct column name
  }
  else {
    # Try a more general lookup, e.g., the third column
    if (ncol(df) > 2 && !paste0("predicted_label_", method_name_pred_col) %in% colnames(df)) {
      message(paste("Assuming third column of", deparse(substitute(df)), "is the prediction. Renaming to 'predicted_label_", method_name_pred_col, "'.", sep=""))
      colnames(df)[3] <- paste0("predicted_label_", method_name_pred_col)
    } else if (!paste0("predicted_label_", method_name_pred_col) %in% colnames(df)){
      stop(paste("Could not find or unambiguously identify prediction column for", method_name_pred_col))
    }
  }
  return(df[, c("cell_id", "true_label", paste0("predicted_label_", method_name_pred_col))])
}

singleR_results_std <- standardize_colnames(singleR_results_raw, "singleR")
scanvi_results_std <- standardize_colnames(scanvi_results_raw, "scanvi")
scmap_results_std <- standardize_colnames(scmap_results_raw, "scmap")


# Find common cell IDs across all methods
common_cell_ids <- Reduce(intersect, list(
  singleR_results_std$cell_id,
  scanvi_results_std$cell_id,
  scmap_results_std$cell_id
))
message(sprintf("Found %d common cell IDs across SingleR, scANVI, and scmap.", length(common_cell_ids)))

if (length(common_cell_ids) == 0) {
  stop("No common cell IDs found across all methods. Check input files and cell ID consistency.")
}

# Merge results, keeping only common cells, and sort by common_cell_ids
# Use a chain of left_join, starting with singleR, then adding scanvi and scmap
# (assuming true_label is consistent across all files for common cells, we take it from singleR)
results_merged <- singleR_results_std %>%
  filter(cell_id %in% common_cell_ids) %>%
  select(cell_id, true_label, predicted_label_singleR) %>%
  arrange(match(cell_id, common_cell_ids))

results_merged <- results_merged %>%
  left_join(
    scanvi_results_std %>%
      filter(cell_id %in% common_cell_ids) %>%
      select(cell_id, predicted_label_scanvi) %>%
      arrange(match(cell_id, common_cell_ids)),
    by = "cell_id"
  )

results_merged <- results_merged %>%
  left_join(
    scmap_results_std %>%
      filter(cell_id %in% common_cell_ids) %>%
      select(cell_id, predicted_label_scmap) %>%
      arrange(match(cell_id, common_cell_ids)),
    by = "cell_id"
  )

# Validate row count after merging
if (nrow(results_merged) != length(common_cell_ids)) {
  stop("Error in merging prediction results. Row count mismatch.")
}
# Validate cell ID order
if (!all(results_merged$cell_id == common_cell_ids)) {
  stop("Cell ID order mismatch after merging.")
}

message("Prediction data successfully loaded, standardized, and merged for common cells.")

true_labels_vec <- results_merged$true_label
singleR_pred_vec <- results_merged$predicted_label_singleR
scanvi_pred_vec <- results_merged$predicted_label_scanvi
scmap_pred_vec <- results_merged$predicted_label_scmap

# Ensure labels are factors with consistent levels
all_unique_labels <- unique(c(
  as.character(true_labels_vec),
  as.character(singleR_pred_vec),
  as.character(scanvi_pred_vec),
  as.character(scmap_pred_vec)
))
# Remove NA (if there are unclassified cells, they will be handled during metric calculation)
all_unique_labels <- sort(na.omit(all_unique_labels))

message(sprintf("Found %d unique cell types (excluding NA): %s", length(all_unique_labels), paste(all_unique_labels, collapse = ", ")))

true_labels_factor <- factor(true_labels_vec, levels = all_unique_labels)
singleR_pred_factor <- factor(singleR_pred_vec, levels = all_unique_labels)
scanvi_pred_factor <- factor(scanvi_pred_vec, levels = all_unique_labels)
scmap_pred_factor <- factor(scmap_pred_vec, levels = all_unique_labels)


# --- 2. Calculate performance metrics ---
message("Calculating performance metrics using MLmetrics and vcd...")

calculate_metrics <- function(pred_labels_factor, true_labels_factor, method_name) {
  # pred_labels_factor and true_labels_factor already use all_unique_labels as levels
  
  # Handle NA predictions: MLmetrics usually doesn't directly handle NA in factors.
  # Here, NA predictions are handled by MLmetrics::Accuracy (treats them as incorrect).
  # For other metrics, we often compare non-NA predictions to non-NA true labels.
  # table() will count NA as an '<NA>' level (if useNA = "ifany" or "always")
  
  # Calculate confusion matrix (using table, which can handle NA in factors if useNA is set)
  cm <- table(Predicted = pred_labels_factor, True = true_labels_factor, useNA = "ifany")
  
  # Overall accuracy (Accuracy from MLmetrics handles NA in y_pred by treating them as incorrect)
  # For this to work best, NA should not be in levels of true_labels_factor unless it's a valid true label
  accuracy <- Accuracy(y_pred = pred_labels_factor, y_true = true_labels_factor)
  
  kappa_val <- NA
  tryCatch({
    # vcd::Kappa expects a confusion matrix without NA row/column names, so we only use non-NA levels for construction
    cm_for_kappa <- table(Predicted = factor(pred_labels_factor, levels=all_unique_labels),
                          True = factor(true_labels_factor, levels=all_unique_labels))
    if (nrow(cm_for_kappa) > 0 && ncol(cm_for_kappa) > 0 && sum(cm_for_kappa) > 0) {
      kappa_result <- vcd::Kappa(cm_for_kappa)
      kappa_val <- kappa_result$Unweighted[1]
    } else {
      warning(paste("Kappa calculation skipped for", method_name, "due to empty confusion matrix after NA removal for Kappa."))
    }
  }, error = function(e) {
    warning(paste("Could not calculate Kappa for", method_name, ":", e$message))
  })
  
  classes <- all_unique_labels # Use predefined non-NA levels
  precision <- numeric(length(classes))
  recall <- numeric(length(classes))
  f1 <- numeric(length(classes))
  names(precision) <- names(recall) <- names(f1) <- classes
  
  for (i in seq_along(classes)) {
    class_val <- classes[i]
    # MLmetrics requires numeric 0/1, the positive parameter specifies which is the positive class
    # For Precision, Recall, F1, we only consider non-NA predictions and true labels
    true_binary <- ifelse(as.character(true_labels_factor) == class_val & !is.na(true_labels_factor), 1, 0)
    pred_binary <- ifelse(as.character(pred_labels_factor) == class_val & !is.na(pred_labels_factor), 1, 0)
    
    # Ensure at least one positive case (true or predicted) to avoid division by zero
    if (sum(true_binary) > 0 || sum(pred_binary > 0)) {
      # If positive = "1" causes issues (e.g., no 1s present), try using positive = 1
      precision[i] <- tryCatch(Precision(y_true = true_binary, y_pred = pred_binary, positive = 1), error = function(e) NA)
      recall[i]    <- tryCatch(Recall(y_true = true_binary, y_pred = pred_binary, positive = 1), error = function(e) NA)
      f1[i]        <- tryCatch(F1_Score(y_true = true_binary, y_pred = pred_binary, positive = 1), error = function(e) NA)
    } else {
      precision[i] <- NA
      recall[i] <- NA
      f1[i] <- NA
    }
  }
  
  # For overall metrics, usually based on non-NA comparisons
  num_correct <- sum(as.character(pred_labels_factor) == as.character(true_labels_factor), na.rm = TRUE)
  num_compared <- sum(!is.na(pred_labels_factor) & !is.na(true_labels_factor)) # Only compare where both are non-NA
  accuracy_manual_no_na <- if(num_compared > 0) num_correct / num_compared else NA
  
  metrics_df <- data.frame(
    Method = method_name,
    Class = classes,
    Precision = precision,
    Recall = recall,
    F1 = f1,
    stringsAsFactors = FALSE
  )
  
  overall_df <- data.frame(
    Method = method_name,
    Accuracy = accuracy, # MLmetrics Accuracy
    Accuracy_NoNA = accuracy_manual_no_na, # Manual calculation, excluding cases where prediction or true label is NA
    Kappa = kappa_val,
    Num_Common_Cells = length(true_labels_factor), # Total number of common cells
    Num_Predicted_Not_NA = sum(!is.na(pred_labels_factor)), # Number of non-NA predictions for this method
    stringsAsFactors = FALSE
  )
  
  # When returning the confusion matrix, also use only non-NA levels
  cm_no_na <- table(Predicted = factor(pred_labels_factor, levels=all_unique_labels),
                    True = factor(true_labels_factor, levels=all_unique_labels))
  
  return(list(confusion_matrix = cm_no_na, metrics_class = metrics_df, metrics_overall = overall_df))
}


singleR_metrics_list <- calculate_metrics(singleR_pred_factor, true_labels_factor, "SingleR")
scanvi_metrics_list <- calculate_metrics(scanvi_pred_factor, true_labels_factor, "scANVI")
scmap_metrics_list <- calculate_metrics(scmap_pred_factor, true_labels_factor, "scmap")

overall_metrics <- rbind(
  singleR_metrics_list$metrics_overall,
  scanvi_metrics_list$metrics_overall,
  scmap_metrics_list$metrics_overall
)

class_metrics <- rbind(
  singleR_metrics_list$metrics_class,
  scanvi_metrics_list$metrics_class,
  scmap_metrics_list$metrics_class
)
class_metrics$F1[is.na(class_metrics$F1)] <- 0 # If F1 is NA (e.g., a class is not present in true and predicted), set to 0

write.csv(overall_metrics, file.path(results_dir, "overall_metrics.csv"), row.names = FALSE)
write.csv(class_metrics, file.path(results_dir, "class_metrics.csv"), row.names = FALSE)
message("Performance metrics saved.")

# --- 3. Visualization: Confusion Matrix Heatmaps ---
plot_confusion_matrix <- function(cm, method_name) {
  cm_df <- as.data.frame(as.table(cm))
  colnames(cm_df) <- c("Predicted", "True", "Count")
  # Ensure factor levels are consistent to prevent ggplot auto-sorting
  cm_df$Predicted <- factor(cm_df$Predicted, levels = rownames(cm))
  cm_df$True <- factor(cm_df$True, levels = colnames(cm))
  
  p <- ggplot(cm_df, aes(x = Predicted, y = True, fill = log1p(Count))) +
    geom_tile(color = "grey50") +
    geom_text(aes(label = Count), color = "black", size = 3) +
    scale_fill_gradient(low = "lightblue", high = "darkblue", name = "log1p(Count)") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(angle = 0, hjust = 1),
          plot.title = element_text(hjust = 0.5)) +
    labs(title = paste(method_name, "Confusion Matrix"), x = "Predicted Label", y = "True Label") +
    coord_fixed()
  return(p)
}

p_singleR_cm <- plot_confusion_matrix(singleR_metrics_list$confusion_matrix, "SingleR")
p_scanvi_cm <- plot_confusion_matrix(scanvi_metrics_list$confusion_matrix, "scANVI")
p_scmap_cm <- plot_confusion_matrix(scmap_metrics_list$confusion_matrix, "scmap")

ggsave(file.path(results_dir, "singleR_confusion_matrix.png"), p_singleR_cm, width = max(7, length(all_unique_labels)*0.8), height = max(6, length(all_unique_labels)*0.7))
ggsave(file.path(results_dir, "scanvi_confusion_matrix.png"), p_scanvi_cm, width = max(7, length(all_unique_labels)*0.8), height = max(6, length(all_unique_labels)*0.7))
ggsave(file.path(results_dir, "scmap_confusion_matrix.png"), p_scmap_cm, width = max(7, length(all_unique_labels)*0.8), height = max(6, length(all_unique_labels)*0.7))
message("Confusion matrices saved.")

# --- 4. Visualization: UMAP Embeddings ---
message("Generating UMAP visualizations...")
sce_query_file <- file.path(results_dir, "sce_query_processed.rds")
if (!file.exists(sce_query_file)) {
  stop(paste("sce_query_processed.rds not found in results_dir:", sce_query_file,
             "\nThis file is needed for UMAP. Ensure it's saved by the R data processing script."))
}
sce_query_original <- readRDS(sce_query_file)

# Subset sce_query_original to match common_cell_ids
# and ensure its column order matches common_cell_ids
cells_in_sce_and_common <- intersect(colnames(sce_query_original), common_cell_ids)
if (length(cells_in_sce_and_common) != length(common_cell_ids)) {
  message(paste("Warning: Not all common_cell_ids (", length(common_cell_ids),
                ") were found in sce_query_processed.rds (", length(cells_in_sce_and_common), " found).",
                "UMAP will be based on this intersection."))
  # If there are indeed ID mismatches, we need to redefine common_cell_ids as this intersection and re-filter results_merged
  common_cell_ids_for_umap <- cells_in_sce_and_common
  results_merged_for_umap <- results_merged %>% filter(cell_id %in% common_cell_ids_for_umap) %>%
    arrange(match(cell_id, common_cell_ids_for_umap))
  if(nrow(results_merged_for_umap) == 0) stop("No cells left for UMAP after intersecting with sce_query_processed.rds.")
} else {
  common_cell_ids_for_umap <- common_cell_ids
  results_merged_for_umap <- results_merged # already sorted by common_cell_ids
}

sce_query_aligned <- sce_query_original[, common_cell_ids_for_umap]

if (!all(colnames(sce_query_aligned) == common_cell_ids_for_umap)) { # Double check order
  sce_query_aligned <- sce_query_aligned[, common_cell_ids_for_umap]
}
message(paste("sce_query object subsetted and aligned to", ncol(sce_query_aligned), "cells for UMAP."))

logcounts_query_aligned <- assay(sce_query_aligned, "logcounts")
if (is.null(logcounts_query_aligned)) stop("'logcounts' not found in aligned query SCE for UMAP.")

set.seed(123)
umap_config <- umap.defaults
# Ensure n_neighbors is valid (less than or equal to number of samples - 1)
# umap input is features x cells, t(logcounts) is cells x features
# so nrow(t(logcounts_query_aligned)) is the number of cells
num_cells_for_umap <- ncol(logcounts_query_aligned)
umap_config$n_neighbors <- min(15, num_cells_for_umap - 1)
umap_config$min_dist <- 0.1

if (umap_config$n_neighbors <= 1 && num_cells_for_umap > 1) { # if only 2 cells, n_neighbors=1
  message("Warning: n_neighbors for UMAP is very small due to low cell count for UMAP. Result might be trivial.")
} else if (num_cells_for_umap <= 1) {
  stop("Cannot run UMAP with <= 1 cell.")
}

message(paste("Running UMAP with n_neighbors =", umap_config$n_neighbors, "and min_dist =", umap_config$min_dist))
umap_result_layout <- umap::umap(t(logcounts_query_aligned), config = umap_config)$layout # umap input is cells x features

# Construct umap_df, all columns based on the order of common_cell_ids_for_umap
umap_df <- data.frame(
  cell_id = common_cell_ids_for_umap,
  UMAP1 = umap_result_layout[, 1],
  UMAP2 = umap_result_layout[, 2],
  True_Label = factor(results_merged_for_umap$true_label, levels = all_unique_labels),
  SingleR_Pred = factor(results_merged_for_umap$predicted_label_singleR, levels = all_unique_labels),
  scANVI_Pred = factor(results_merged_for_umap$predicted_label_scanvi, levels = all_unique_labels),
  scmap_Pred = factor(results_merged_for_umap$predicted_label_scmap, levels = all_unique_labels)
)

plot_umap <- function(df, label_col, title_str) {
  p <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = .data[[label_col]])) +
    geom_point(size = 0.5, alpha = 0.7) +
    theme_minimal(base_size = 10) +
    guides(color = guide_legend(override.aes = list(size=3), title = "Cell Label")) +
    labs(title = title_str) +
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5))
  return(p)
}

p_umap_true <- plot_umap(umap_df, "True_Label", "True Label UMAP")
p_umap_singleR <- plot_umap(umap_df, "SingleR_Pred", "SingleR Prediction UMAP")
p_umap_scanvi <- plot_umap(umap_df, "scANVI_Pred", "scANVI Prediction UMAP")
p_umap_scmap <- plot_umap(umap_df, "scmap_Pred", "scmap Prediction UMAP")

p_umap_combined <- (p_umap_true | p_umap_singleR) / (p_umap_scanvi | p_umap_scmap) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave(file.path(results_dir, "umap_predictions_combined.png"), p_umap_combined, width = 12, height = 11)
message("UMAP visualizations saved.")

# --- 5. Visualization: Performance Metric Bar Plots ---
message("Generating performance bar plots...")
# Use Accuracy (MLmetrics version, handles NA predictions as incorrect)
p_accuracy <- ggplot(overall_metrics, aes(x = reorder(Method, -Accuracy), y = Accuracy, fill = Method)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  geom_text(aes(label=sprintf("%.3f", Accuracy)), vjust=-0.3, size=3.5) +
  theme_minimal(base_size = 10) +
  labs(title = "Overall Accuracy Comparison (NA predictions counted as incorrect)", y = "Accuracy", x = "Method") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))

p_kappa <- ggplot(overall_metrics, aes(x = reorder(Method, -Kappa), y = Kappa, fill = Method)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  geom_text(aes(label=sprintf("%.3f", Kappa)), vjust=-0.3, size=3.5) +
  theme_minimal(base_size = 10) +
  labs(title = "Kappa Coefficient Comparison (based on non-NA predictions)", y = "Kappa Coefficient", x = "Method") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))

p_f1 <- ggplot(class_metrics %>% filter(!is.na(F1)), aes(x = Class, y = F1, fill = Method)) + # Filter out classes where F1 is NA
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top") +
  labs(title = "F1 Score Comparison by Class", y = "F1 Score", x = "Cell Type") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))

num_classes_plot <- length(unique(class_metrics$Class[!is.na(class_metrics$F1)]))
f1_plot_width <- max(8, num_classes_plot * 0.5 + 3)

ggsave(file.path(results_dir, "accuracy_comparison.png"), p_accuracy, width = 6, height = 5)
ggsave(file.path(results_dir, "kappa_comparison.png"), p_kappa, width = 6, height = 5)
ggsave(file.path(results_dir, "f1_comparison.png"), p_f1, width = f1_plot_width, height = 6)
message("Performance bar plots saved.")


# --- 6. Output Detailed Text File ---
message("Generating detailed benchmark report...")
report_file <- file.path(results_dir, "benchmark_report.txt")
cat("Benchmark Analysis Report (Detailed)\n", file = report_file)
cat(paste("Report generated on:", Sys.time(), "\n"), file = report_file, append = TRUE)
cat(paste("Number of common cells analyzed (across all methods):", length(common_cell_ids), "\n"), file = report_file, append = TRUE)
cat(paste("Number of cells used for UMAP (intersection with sce_query_processed.rds):", length(common_cell_ids_for_umap), "\n"), file = report_file, append = TRUE)
cat(paste("Cell types considered (excluding NA):", paste(all_unique_labels, collapse = ", "), "\n\n"), file = report_file, append = TRUE)

cat("--- Overall Performance Metrics ---\n", file = report_file, append = TRUE)
write.table(overall_metrics, report_file, append = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
cat("\n\n", file = report_file, append = TRUE)

cat("--- Class-wise Performance Metrics (Precision, Recall, F1) ---\n", file = report_file, append = TRUE)
class_metrics_wide <- class_metrics %>%
  tidyr::pivot_wider(names_from = Method, values_from = c(Precision, Recall, F1), names_vary = "slowest")
write.table(class_metrics_wide, report_file, append = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
cat("\n\n", file = report_file, append = TRUE)

cat("--- Confusion Matrices (Predicted vs True; based on non-NA labels) ---\n", file = report_file, append = TRUE)
cat("\nSingleR Confusion Matrix:\n", file = report_file, append = TRUE)
write.table(as.matrix(singleR_metrics_list$confusion_matrix), report_file, append = TRUE, row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
cat("\n\n", file = report_file, append = TRUE)

cat("scANVI Confusion Matrix:\n", file = report_file, append = TRUE)
write.table(as.matrix(scanvi_metrics_list$confusion_matrix), report_file, append = TRUE, row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
cat("\n\n", file = report_file, append = TRUE)

cat("scmap Confusion Matrix:\n", file = report_file, append = TRUE)
write.table(as.matrix(scmap_metrics_list$confusion_matrix), report_file, append = TRUE, row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
cat("\n", file = report_file, append = TRUE)

message("Benchmark analysis complete. Results and visualizations saved in ", results_dir)
message("Detailed text report saved to: ", report_file)
