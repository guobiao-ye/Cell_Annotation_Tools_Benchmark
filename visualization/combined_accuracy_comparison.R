# Load necessary libraries
library(ggplot2)
library(dplyr)

# --- Accuracy data for Combined dataset ---
# Format: Proportion (training set proportion), Method, Accuracy
combined_data <- data.frame(
  Proportion_Text = factor(c(
    "0.75", "0.75", "0.75",
    "0.25", "0.25", "0.25",
    "0.125", "0.125", "0.125",
    "0.025", "0.025", "0.025"
  ), levels = c("0.025", "0.125", "0.25", "0.75")), # For ordered X-axis in plotting
  Method = factor(rep(c("SingleR", "scANVI", "scmap"), times = 4),
                  levels = c("SingleR", "scANVI", "scmap")),
  Accuracy = c(
    # 0.75 Training (Combined_0.75_training_benchmark_report.txt)
    0.978062157221207, 1, 0.99908592321755,
    # 0.25 Training (Combined_0.25_training_benchmark_report.txt)
    0.981707317073171, 0.99969512195122, 0.985365853658537,
    # 0.125 Training (Combined_0.125_training_benchmark_report.txt)
    0.981970211654037, 0.999738698719624, 0.974653775803501,
    # 0.025 Training (Combined_0.025_training_benchmark_report.txt)
    0.956848030018762, 0.980300187617261, 0.895637898686679
  )
)

# Convert Proportion_Text to numeric for plotting, but keep factor for ordering
combined_data$Proportion_Numeric <- as.numeric(as.character(combined_data$Proportion_Text))

# Plot bar chart
p_combined <- ggplot(combined_data, aes(x = Proportion_Text, y = Accuracy, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
  geom_text(aes(label = sprintf("%.3f", Accuracy)),
            position = position_dodge(width = 0.9),
            vjust = -0.25, size = 3) +
  scale_fill_manual(values = c("SingleR" = "#1f77b4", "scANVI" = "#ff7f0e", "scmap" = "#2ca02c")) +
  labs(title = "Accuracy Comparison on Combined (10X + CEL-Seq2) Dataset",
       x = "Reference Set Proportion",
       y = "Overall Accuracy",
       fill = "Method") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "top"
  ) +
  ylim(0, 1.05) # Ensure labels do not exceed bounds

# Display the plot
print(p_combined)

# Save the plot
ggsave("combined_accuracy_comparison.png", plot = p_combined, width = 10, height = 6, dpi = 300)