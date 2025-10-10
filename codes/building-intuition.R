library(dtw)
library(ggplot2)
library(reshape2)
library(gridExtra)
# this file shows examples and DTW alignments (as in Figure 2 in the preprint). 
# Generate synthetic time series data
set.seed(123)
time_series_1 <- sin(seq(0, 2 * pi, length.out = 100)) + rnorm(100, sd = 0.1)
time_series_2 <- sin(seq(0, 2 * pi, length.out = 120)) + rnorm(120, sd = 0.1)
time_series_3 <- cos(seq(0, 2 * pi, length.out = 100)) + rnorm(100, sd = 0.1)
time_series_4 <- sin(seq(0, 2 * pi, length.out = 100)) + 0.3 + rnorm(100, sd = 0.1)
time_series_shifted <- c(rep(NA, 20), time_series_1[1:(length(time_series_1) - 20)])

# Remove NA values from the shifted series for DTW computation
time_series_shifted <- na.omit(time_series_shifted)

# Function to compute DTW and prepare data for plotting
compute_dtw <- function(ts1, ts2, series1_name, series2_name) {
  alignment <- dtw(ts1, ts2, keep = TRUE)
  warp_path <- data.frame(index1 = alignment$index1, index2 = alignment$index2)
  
  ts1_df <- data.frame(index = 1:length(ts1), value = ts1, series = series1_name)
  ts2_df <- data.frame(index = 1:length(ts2), value = ts2, series = series2_name)
  
  combined_df <- rbind(ts1_df, ts2_df)
  list(combined_df = combined_df, warp_path = warp_path, ts1 = ts1, ts2 = ts2)
}

# Compute DTW alignments
dtw_results_1_2 <- compute_dtw(time_series_1, time_series_2, "Series 1", "Series 2")
dtw_results_1_3 <- compute_dtw(time_series_1, time_series_3, "Series 1", "Series 3")
dtw_results_1_4 <- compute_dtw(time_series_1, time_series_4, "Series 1", "Series 4")
dtw_results_1_shifted <- compute_dtw(time_series_1, time_series_shifted, "Series 1", "Series 1 Shifted")

# Plot time series and DTW alignments
plot_dtw <- function(dtw_results, title, show_y_label) {
  p1 <- ggplot() +
    geom_line(data = dtw_results$combined_df, aes(x = index, y = value, color = series)) +
    geom_segment(data = dtw_results$warp_path, aes(x = index1, xend = index2, y = dtw_results$ts1[index1], yend = dtw_results$ts2[index2]), alpha = 0.2, linetype = "dotted") +
    labs(title = title, x = "Time Index", y = "Value", color = " ") +
    theme_minimal()
  
  p2 <- ggplot(dtw_results$warp_path, aes(x = index1, y = index2)) +
    geom_line() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    labs(title = " ", x = "Index in Series 1", y = if (show_y_label) "Index in Comparison Series" else NULL) +
    theme_minimal()
  
  list(p1, p2)
}

# Plot all DTW alignments and manage labels
plot1 <- plot_dtw(dtw_results_1_2, " ", FALSE)
plot2 <- plot_dtw(dtw_results_1_3, " ", TRUE)
plot3 <- plot_dtw(dtw_results_1_4, " ", FALSE)
plot4 <- plot_dtw(dtw_results_1_shifted, " ", FALSE)

# Save the arranged plots as a PDF
pdf("DTW_Alignment_PlotsXXX.pdf", width = 8, height = 8)

grid.arrange(plot1[[1]], plot1[[2]], 
             plot2[[1]], plot2[[2]], 
             plot3[[1]], plot3[[2]], 
             plot4[[1]], plot4[[2]], 
             ncol = 2)

dev.off()
