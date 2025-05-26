# Load necessary libraries
library(ggplot2)
library(dplyr)
library(scales) # For percentage labels

# --- Execution time data for Combined dataset (seconds) ---
combined_time_data_ext <- data.frame(
  Proportion_Numeric = c(
    0.75, 0.75, 0.75,        
    0.50, 0.50, 0.50,        
    0.25, 0.25, 0.25,        
    0.125, 0.125, 0.125,     
    0.05, 0.05, 0.05,        
    0.025, 0.025, 0.025,     
    0.0125, 0.0125, 0.0125   
  ),
  Method = factor(rep(c("SingleR", "scANVI", "scmap"), times = 7),
                  levels = c("SingleR", "scANVI", "scmap")),
  Total_Execution_Time = c(
    # 0.75
    7.482671,  # SingleR
    (39.095456 + 5.177625 + 0.032472), # scANVI
    3.710033,   # scmap
    
    # 0.50
    7.138324,       # SingleR
    (26.246883 + 3.364682 + 0.048358), # scANVI
    4.853243,       # scmap
    
    # 0.25
    6.525334,  # SingleR
    (14.530908 + 2.955436 + 0.065918), # scANVI
    6.011979,   # scmap
    
    # 0.125
    3.863534,  # SingleR
    (8.244467 + 1.864307 + 1.206973), # scANVI
    6.550142,   # scmap
    
    # 0.05
    2.462892,       # SingleR
    (5.173529 + 1.482316 + 0.473456),  # scANVI
    7.948424,       # scmap
    
    # 0.025
    1.938152,  # SingleR
    (3.807783 + 1.840196 + 0.080123), # scANVI
    10.037626,  # scmap
    
    # 0.0125
    1.442472,       # SingleR
    (2.428911 + 1.215779 + 0.187492), # scANVI
    11.063792       # scmap
  )
)

# Plot line chart
p_time_combined_ext <- ggplot(combined_time_data_ext,
                              aes(x = Proportion_Numeric, y = Total_Execution_Time, group = Method, color = Method)) +
  geom_line(linewidth = 1.2) + # Thicker lines
  geom_point(size = 3.5, aes(shape = Method)) + # Add points with different shapes
  scale_color_manual(values = c("SingleR" = "#1f77b4", "scANVI" = "#ff7f0e", "scmap" = "#2ca02c")) +
  scale_shape_manual(values = c("SingleR" = 16, "scANVI" = 17, "scmap" = 15)) + # Circle, triangle, square
  labs(title = "Execution Time Comparison on Combined (10X + CEL-Seq2) Dataset",
       x = "Reference Set Proportion",
       y = "Total Execution Time (seconds)",
       color = "Method",
       shape = "Method") +
  theme_minimal(base_size = 15) + # Increase base font size
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    legend.position = "top",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.major = element_line(colour = "grey85"), # Major grid lines
    panel.grid.minor = element_line(colour = "grey90")  # Minor grid lines
  ) +
  scale_x_continuous(
    breaks = seq(0, 0.80, by = 0.10), # Regular intervals, e.g., every 10%
    labels = scales::percent_format(accuracy = 1), # X-axis in percentage
    limits = c(0, 0.80) # Set X-axis range to ensure all points are visible
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x), # Standard log scale
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  annotation_logticks(sides = "l") # Add log ticks on the left Y-axis

# Display the plot
print(p_time_combined_ext)

# Save the plot if needed
ggsave("combined_time_comparison.png", plot = p_time_combined_ext, width = 12, height = 7, dpi = 300)
