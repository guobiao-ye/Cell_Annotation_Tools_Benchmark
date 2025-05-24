# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr) # For pivot_longer and other data tidying
library(scales) # For percentage and scientific notation labels

# --- Estimated memory/GPU memory data for Combined dataset (GB) ---
# 'ResourceType' column distinguishes CPU_RAM and GPU_VRAM
# GPU_VRAM is only filled for scANVI, others are NA
combined_resource_data <- data.frame(
  Proportion_Numeric = c(
    # 75%
    0.75, 0.75, 0.75, 0.75, 0.75,
    # 50%
    0.50, 0.50, 0.50, 0.50, 0.50,
    # 25%
    0.25, 0.25, 0.25, 0.25, 0.25,
    # 12.5%
    0.125, 0.125, 0.125, 0.125, 0.125,
    # 5%
    0.05, 0.05, 0.05, 0.05, 0.05,
    # 2.5%
    0.025, 0.025, 0.025, 0.025, 0.025,
    # 1.25%
    0.0125, 0.0125, 0.0125, 0.0125, 0.0125
  ),
  Method = factor(rep(c("SingleR", "scANVI", "scANVI", "scmap", "scmap"), times = 7), # scANVI repeated for RAM and VRAM
                  levels = c("SingleR", "scANVI", "scmap")),
  ResourceType = factor(rep(c("CPU_RAM", "CPU_RAM", "GPU_VRAM", "CPU_RAM", "CPU_RAM_placeholder_scmap_vram"), times = 7), # Placeholder for scmap
                        levels = c("CPU_RAM", "GPU_VRAM", "CPU_RAM_placeholder_scmap_vram")),
  Value = c(
    # 0.75
    3.15, 12.55, 6.75, 2.95, NA, # SingleR CPU, scANVI CPU, scANVI GPU, scmap CPU, scmap GPU (NA)
    # 0.50
    3.05, 9.35, 4.65, 2.75, NA,
    # 0.25
    2.25, 6.15, 3.25, 2.10, NA,
    # 0.125
    1.55, 4.20, 2.15, 1.65, NA,
    # 0.05
    1.10, 2.75, 1.45, 1.40, NA,
    # 0.025
    0.85, 2.15, 1.10, 1.25, NA,
    # 0.0125
    0.70, 1.50, 0.85, 1.30, NA
  )
)

# Remove scmap's GPU placeholder rows (since scmap does not use GPU)
combined_resource_data <- combined_resource_data %>%
  filter(!(Method == "scmap" & ResourceType == "CPU_RAM_placeholder_scmap_vram"))
# Update ResourceType factor levels
combined_resource_data$ResourceType <- factor(combined_resource_data$ResourceType, levels = c("CPU_RAM", "GPU_VRAM"))

# 1. Line plot visualization
p_resource_line_combined <- ggplot(combined_resource_data,
                                   aes(x = Proportion_Numeric, y = Value,
                                       group = interaction(Method, ResourceType), # Group by method and resource type
                                       color = Method, linetype = ResourceType)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.5, aes(shape = Method)) +
  scale_color_manual(values = c("SingleR" = "#1f77b4", "scANVI" = "#ff7f0e", "scmap" = "#2ca02c"), name = "Method") +
  scale_linetype_manual(values = c("CPU_RAM" = "solid", "GPU_VRAM" = "dashed"), name = "Resource Type") + # Solid for CPU, dashed for GPU
  scale_shape_manual(values = c("SingleR" = 16, "scANVI" = 17, "scmap" = 15), name = "Method") +
  labs(title = "Estimated Resource Consumption on Combined Dataset",
       x = "Reference Set Proportion (Training)",
       y = "Estimated Peak Consumption (GB, log scale)") +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    legend.position = "top",
    legend.box = "vertical", # Arrange legend vertically
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.major = element_line(colour = "grey85"),
    panel.grid.minor = element_line(colour = "grey90")
  ) +
  scale_x_continuous(
    breaks = unique(combined_resource_data$Proportion_Numeric), # Set breaks based on data points
    labels = scales::percent_format(accuracy = 1),
    limits = c(min(combined_resource_data$Proportion_Numeric)-0.01, max(combined_resource_data$Proportion_Numeric)+0.01) # Adjust X-axis range
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x, n=4), # Automatically determine log scale ticks
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  annotation_logticks(sides = "l")

print(p_resource_line_combined)
# ggsave("combined_resource_line_plot.png", plot = p_resource_line_combined, width = 12, height = 7, dpi = 300)

# 2. Grouped stacked bar plot (scANVI's CPU RAM and GPU VRAM will be stacked, SingleR and scmap only have CPU RAM)
# To stack, ensure scANVI's CPU_RAM and GPU_VRAM are different 'Stage' or 'Component' under the same Method
# Here we only visualize scANVI's stacking, others show their CPU RAM

# Prepare data for stacked plot
# SingleR and scmap only have CPU_RAM, scANVI has CPU_RAM and GPU_VRAM
plot_data_stacked <- combined_resource_data %>%
  mutate(Proportion_Factor = factor(Proportion_Numeric))

# Use a striking color for scANVI's GPU_VRAM
resource_colors_stacked <- c("CPU_RAM" = "lightblue", "GPU_VRAM" = "#e31a1c") # Red for GPU VRAM

p_resource_stacked_combined <- ggplot() +
  # Plot SingleR's CPU RAM
  geom_bar(data = filter(plot_data_stacked, Method == "SingleR", ResourceType == "CPU_RAM"),
           aes(x = Proportion_Factor, y = Value, fill = ResourceType),
           stat = "identity", position = "stack", width = 0.25) +
  # Plot scmap's CPU RAM
  geom_bar(data = filter(plot_data_stacked, Method == "scmap", ResourceType == "CPU_RAM"),
           aes(x = Proportion_Factor, y = Value, fill = ResourceType),
           stat = "identity", position = "stack", width = 0.25) +
  # Plot scANVI's CPU RAM and GPU VRAM (stacked)
  geom_bar(data = filter(plot_data_stacked, Method == "scANVI"),
           aes(x = Proportion_Factor, y = Value, fill = ResourceType),
           stat = "identity", position = "stack", width = 0.25) +
  facet_wrap(~ Method, scales = "free_x", ncol = 3) +
  scale_fill_manual(values = resource_colors_stacked, name = "Resource Type") +
  labs(title = "Estimated Resource Consumption on Combined (10X + CEL-Seq2) Dataset",
       x = "Reference Set Proportion (Training)",
       y = "Estimated Peak Consumption (GB, log scale)") +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_text(size = 14),
    legend.position = "top",
    strip.text = element_text(face = "bold", size = 12),
    panel.spacing = unit(1, "lines")
  ) +
  scale_x_discrete(labels = function(x) paste0(as.numeric(x) * 100, "%")) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  annotation_logticks(sides = "l")

print(p_resource_stacked_combined)
ggsave("combined_resource_stacked_plot.png", plot = p_resource_stacked_combined, width = 14, height = 7, dpi = 300)