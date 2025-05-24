# Load necessary libraries
library(ggplot2)
library(dplyr)

# --- Accuracy data for CelSeq2_5cl dataset ---
# Format: Proportion (training set proportion), Method, Accuracy
celseq2_data <- data.frame(
  Proportion_Text = factor(c(
    "0.25", "0.25", "0.25",
    "0.125", "0.125", "0.125",
    "0.05", "0.05", "0.05",
    "0.025", "0.025", "0.025"
  ), levels = c("0.025", "0.05", "0.125", "0.25")), # For ordered X-axis in plotting
  Method = factor(rep(c("SingleR", "scANVI", "scmap"), times = 4),
                  levels = c("SingleR", "scANVI", "scmap")),
  Accuracy = c(
    # 0.25 Training (CelSeq2_5cl_0.25_training_benchmark_report.txt)
    0.936915887850467, 1, 0.883177570093458,
    # 0.125 Training (CelSeq2_5cl_0.125_training_benchmark_report.txt)
    0.925851703406814, 0.997995991983968, 0.803607214428858,
    # 0.05 Training (CelSeq2_5cl_0.05_training_benchmark_report.txt)
    0.994464944649446, 0.623616236162362, 0.601476014760148,
    # 0.025 Training (CelSeq2_5cl_0.025_training_benchmark_report.txt)
    0.982014388489209, 0.811151079136691, 0.309352517985612
  )
)

# Convert Proportion_Text to numeric for plotting, but keep factor for ordering
celseq2_data$Proportion_Numeric <- as.numeric(as.character(celseq2_data$Proportion_Text))

# Plot bar chart
p_celseq2 <- ggplot(celseq2_data, aes(x = Proportion_Text, y = Accuracy, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
  geom_text(aes(label = sprintf("%.3f", Accuracy)),
            position = position_dodge(width = 0.9),
            vjust = -0.25, size = 3) +
  scale_fill_manual(values = c("SingleR" = "#1f77b4", "scANVI" = "#ff7f0e", "scmap" = "#2ca02c")) +
  labs(title = "Accuracy Comparison on CEL-Seq2 Dataset",
       x = "Reference Set Proportion",
       y = "Overall Accuracy",
       fill = "Method") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5), # Change to angle = 45, hjust = 1 if labels overlap
    legend.position = "top"
  ) +
  ylim(0, 1.05) # Ensure labels do not exceed bounds

# Display the plot
print(p_celseq2)

# Save the plot
ggsave("celseq2_accuracy_comparison.png", plot = p_celseq2, width = 10, height = 6, dpi = 300)