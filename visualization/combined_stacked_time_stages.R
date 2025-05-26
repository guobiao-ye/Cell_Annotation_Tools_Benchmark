# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr) # For pivot_longer
library(scales) # For percentage labels

# --- Execution time data for each stage of Combined dataset (seconds) ---
# Stage: 'Annotation' (for SingleR),
#        'scVI_Train', 'scANVI_Finetune', 'Predict' (for scANVI)
#        'FeatSel', 'Indexing', 'Predict' (for scmap)

combined_stages_time_data <- data.frame(
  Proportion_Numeric = numeric(),
  Method = factor(),
  Stage = factor(),
  Time = numeric()
)

# Helper function to add data
add_method_data <- function(df, prop, method_name, times_list, stage_names_list) {
  for (i in 1:length(stage_names_list)) {
    df <- rbind(df, data.frame(
      Proportion_Numeric = prop,
      Method = method_name,
      Stage = stage_names_list[i],
      Time = times_list[i]
    ))
  }
  return(df)
}

# --- 0.75 Training ---
# SingleR
combined_stages_time_data <- add_method_data(combined_stages_time_data, 0.75, "SingleR",
                                             c(7.482671), c("Annotation"))
# scANVI
combined_stages_time_data <- add_method_data(combined_stages_time_data, 0.75, "scANVI",
                                             c(39.095456, 5.177625, 0.032472),
                                             c("scVI_Train", "scANVI_Finetune", "Predict"))
# scmap
combined_stages_time_data <- add_method_data(combined_stages_time_data, 0.75, "scmap",
                                             c(1.900760, 0.751975, 1.057298),
                                             c("FeatSel", "Indexing", "Predict"))

# --- 0.50 Training ---
# SingleR
combined_stages_time_data <- add_method_data(combined_stages_time_data, 0.50, "SingleR",
                                             c(7.138324), c("Annotation"))
# scANVI
combined_stages_time_data <- add_method_data(combined_stages_time_data, 0.50, "scANVI",
                                             c(26.246883, 3.364682, 0.048358),
                                             c("scVI_Train", "scANVI_Finetune", "Predict"))
# scmap
combined_stages_time_data <- add_method_data(combined_stages_time_data, 0.50, "scmap",
                                             c(1.102461, 0.481881, 3.268901),
                                             c("FeatSel", "Indexing", "Predict"))

# --- 0.25 Training ---
# SingleR
combined_stages_time_data <- add_method_data(combined_stages_time_data, 0.25, "SingleR",
                                             c(6.525334), c("Annotation"))
# scANVI
combined_stages_time_data <- add_method_data(combined_stages_time_data, 0.25, "scANVI",
                                             c(14.530908, 2.955436, 0.065918),
                                             c("scVI_Train", "scANVI_Finetune", "Predict"))
# scmap
combined_stages_time_data <- add_method_data(combined_stages_time_data, 0.25, "scmap",
                                             c(0.268833, 0.255441, 5.487705),
                                             c("FeatSel", "Indexing", "Predict"))

# --- 0.125 Training ---
# SingleR
combined_stages_time_data <- add_method_data(combined_stages_time_data, 0.125, "SingleR",
                                             c(3.863534), c("Annotation"))
# scANVI
combined_stages_time_data <- add_method_data(combined_stages_time_data, 0.125, "scANVI",
                                             c(8.244467, 1.864307, 1.206973),
                                             c("scVI_Train", "scANVI_Finetune", "Predict"))
# scmap
combined_stages_time_data <- add_method_data(combined_stages_time_data, 0.125, "scmap",
                                             c(0.136116, 0.151635, 6.262391),
                                             c("FeatSel", "Indexing", "Predict"))

# --- 0.05 Training ---
# SingleR
combined_stages_time_data <- add_method_data(combined_stages_time_data, 0.05, "SingleR",
                                             c(2.462892), c("Annotation"))
# scANVI
combined_stages_time_data <- add_method_data(combined_stages_time_data, 0.05, "scANVI",
                                             c(5.173529, 1.482316, 0.473456),
                                             c("scVI_Train", "scANVI_Finetune", "Predict"))
# scmap
combined_stages_time_data <- add_method_data(combined_stages_time_data, 0.05, "scmap",
                                             c(0.080242, 0.069825, 7.798357),
                                             c("FeatSel", "Indexing", "Predict"))

# --- 0.025 Training ---
# SingleR
combined_stages_time_data <- add_method_data(combined_stages_time_data, 0.025, "SingleR",
                                             c(1.938152), c("Annotation"))
# scANVI
combined_stages_time_data <- add_method_data(combined_stages_time_data, 0.025, "scANVI",
                                             c(3.807783, 1.840196, 0.080123),
                                             c("scVI_Train", "scANVI_Finetune", "Predict"))
# scmap
combined_stages_time_data <- add_method_data(combined_stages_time_data, 0.025, "scmap",
                                             c(0.124933, 0.113169, 9.799524),
                                             c("FeatSel", "Indexing", "Predict"))

# --- 0.0125 Training ---
# SingleR
combined_stages_time_data <- add_method_data(combined_stages_time_data, 0.0125, "SingleR",
                                             c(1.442472), c("Annotation"))
# scANVI
combined_stages_time_data <- add_method_data(combined_stages_time_data, 0.0125, "scANVI",
                                             c(2.428911, 1.215779, 0.187492),
                                             c("scVI_Train", "scANVI_Finetune", "Predict"))
# scmap
combined_stages_time_data <- add_method_data(combined_stages_time_data, 0.0125, "scmap",
                                             c(0.049246, 0.048892, 10.965654),
                                             c("FeatSel", "Indexing", "Predict"))

# Ensure factor order
combined_stages_time_data$Method <- factor(combined_stages_time_data$Method, levels = c("SingleR", "scmap", "scANVI"))
combined_stages_time_data$Stage <- factor(combined_stages_time_data$Stage, levels = c(
  "Annotation", # SingleR
  "FeatSel", "Indexing", "Predict", # scmap (Predict for scmap comes after Indexing)
  "scVI_Train", "scANVI_Finetune", "Predict_scANVI" # scANVI (renaming scANVI's predict to avoid clash)
))
# Rename scANVI's Predict stage to avoid confusion with scmap's Predict in the legend, if colors are assigned by Stage
combined_stages_time_data$Stage[combined_stages_time_data$Method == "scANVI" & combined_stages_time_data$Stage == "Predict"] <- "Predict_scANVI"
# Update factor levels to include renamed stage
combined_stages_time_data$Stage <- factor(combined_stages_time_data$Stage, levels = c(
  "Annotation", "FeatSel", "Indexing", "Predict",
  "scVI_Train", "scANVI_Finetune", "Predict_scANVI"
))

# Define stage colors (ensure correspondence with Stage factor levels above)
stage_colors <- c(
  "Annotation" = "grey50",      # SingleR
  "FeatSel" = "#fdc086",      # scmap - Feature Selection
  "Indexing" = "#beaed4",     # scmap - Indexing
  "Predict" = "#7fc97f",      # scmap - Prediction
  "scVI_Train" = "#ffff99",   # scANVI - scVI Training (light yellow)
  "scANVI_Finetune" = "#386cb0", # scANVI - Finetuning (deep blue)
  "Predict_scANVI" = "#f0027f"  # scANVI - Prediction (pink)
)

# Plot grouped stacked bar chart
p_stacked_time_combined <- ggplot(combined_stages_time_data,
                                  aes(x = factor(Proportion_Numeric), y = Time, fill = Stage)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  facet_wrap(~ Method, scales = "free_x", ncol = 3) + # One panel per method, X-axis independent to show all proportions
  scale_fill_manual(values = stage_colors, name = "Execution Stage") +
  labs(title = "Execution Time Breakdown on Combined (10X + CEL-Seq2) Dataset by Reference Proportion",
       x = "Reference Set Proportion",
       y = "Execution Time (seconds, log scale)") +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_text(size = 14),
    legend.position = "top",
    strip.text = element_text(face = "bold", size = 12), # Panel titles
    panel.spacing = unit(1, "lines")
  ) +
  scale_x_discrete(labels = function(x) paste0(as.numeric(x) * 100, "%")) + # X-axis in percentage
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  annotation_logticks(sides = "l")

# Display the plot
print(p_stacked_time_combined)

# Save the plot if needed
ggsave("combined_stacked_time_stages.png", plot = p_stacked_time_combined, width = 14, height = 7, dpi = 300)
