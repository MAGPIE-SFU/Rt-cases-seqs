library(perm)
library(dplyr)
library(mgcv)
library(scales)  
#install.packages("dtw")
library(patchwork)
library(philentropy)
library(dtw)
library(ggplot2)
library(gridExtra)

library(magick)
library(cowplot)
library(ggplot2)
library(grid)

# this file reads csv files corresponding to the simulated outbreaks
# and generates plots 


# Function to read CSV files from specified directories and
# match them based on the common number in their names
read_and_match_files <- function(outbreak_dir, trueRt_dir, seqEstim_dir, epiEstim_dir) {
  outbreak_files <- list.files(outbreak_dir, pattern = "\\.csv$", full.names = TRUE)
  outbreak_files <- outbreak_files[!grepl("gt\\.csv$", outbreak_files)] # Exclude files ending with gt.csv
  
  trueRt_files <- list.files(trueRt_dir, pattern = "\\.csv$", full.names = TRUE)
  seqEstim_files <- list.files(seqEstim_dir, pattern = "\\.csv$", full.names = TRUE)
  epiEstim_files <- list.files(epiEstim_dir, pattern = "\\.csv$", full.names = TRUE)

  # Check if any directory is empty
  if (length(outbreak_files) == 0 || length(trueRt_files) == 0 || 
      length(seqEstim_files) == 0 || length(epiEstim_files) == 0) {
    stop("One or more directories are empty. Check file paths and contents.")
  }
  
  # Extract the file numbers (001 to 020) from the file names
  extract_number <- function(filename) {
    num <- sub(".*?(\\d{3}).*", "\\1", basename(filename))
    if (num == filename) return(NA)  # If no number found, return NA
    return(num)
  }
  
  # Group files by their numbers
  grouped_files <- lapply(list(outbreak_files, trueRt_files, seqEstim_files, epiEstim_files), function(file_list) {
    if (length(file_list) == 0) return(NULL)  # Handle empty lists
    extracted_numbers <- sapply(file_list, extract_number)
    extracted_numbers <- na.omit(extracted_numbers)  # Remove NAs
    split(file_list, extracted_numbers)
  })
  
  # Find common keys (file numbers) across all groups
  common_keys <- Reduce(intersect, lapply(grouped_files, names))
  
  if (length(common_keys) == 0) {
    stop("No matching file numbers found across directories. Check file naming consistency.")
  }
  
  # Filter groups that have matching files across all categories
  matched_groups <- lapply(common_keys, function(key) {
    list(outbreak = grouped_files[[1]][[key]][1],
         trueRt = grouped_files[[2]][[key]][1],
         seqEstim = grouped_files[[3]][[key]][1],
         epiEstim = grouped_files[[4]][[key]][1],
         number = key)
  })
  
  return(matched_groups)
}

# Function to process each group of matched files
process_group <- function(group, results_list) {
  # Extract the scenario number
  scenario_number <- group$number
  
  # Print file paths for debugging
  print(paste("Processing files:", group$outbreak, group$trueRt, group$seqEstim, group$epiEstim))
  
  # Check if files exist before reading
  if (!file.exists(group$outbreak)) stop(paste("File not found:", group$outbreak))
  if (!file.exists(group$trueRt)) stop(paste("File not found:", group$trueRt))
  if (!file.exists(group$seqEstim)) stop(paste("File not found:", group$seqEstim))
  if (!file.exists(group$epiEstim)) stop(paste("File not found:", group$epiEstim))
  
  # Read the data from the files
  outbreak_data <- read.csv(group$outbreak)
  true_data <- read.csv(group$trueRt)
  seq_data <- read.csv(group$seqEstim)
  epi_data <- read.csv(group$epiEstim)
  
  # Print column names for debugging
  print("Column names in seq_data:")
  print(colnames(seq_data))
  print("Column names in epi_data:")
  print(colnames(epi_data))
  print("Column names in true_data:")
  print(colnames(true_data))
  print("Column names in outbreak_data:")
  print(colnames(outbreak_data))
  
  # Ensure necessary columns are present
  if (!all(c("med", "lower", "upper", "time") %in% colnames(seq_data))) stop("Missing columns in seq_data")
  if (!all(c("t_end", "Quantile.0.05.R.", "Quantile.0.95.R.", "Median.R.") %in% colnames(epi_data))) stop("Missing columns in epi_data")
  if (!all(c("smoothRt_day_avg", "infectious_time") %in% colnames(true_data))) stop("Missing columns in true_data")
  if (!"incidence" %in% colnames(outbreak_data)) stop("Missing 'incidence' column in outbreak_data")
  
  
  #Truncate true_data and epi_data to match the length of seq_data
  #min_length <- min(nrow(seq_data), nrow(true_data), nrow(epi_data))
  true_data <- true_data#[1:min_length, ]
  epi_data <- epi_data#[1:min_length, ]
  seq_data <- seq_data#[1:min_length, ]
  
  # Apply the warping function
  alignment3 <- dtw(seq_data$med, true_data$smoothRt_day_avg, keep = TRUE)
  alignment4 <- dtw(epi_data$`Median.R.`, true_data$smoothRt_day_avg, keep = TRUE)
  alignment_perfect <- dtw(true_data$smoothRt_day_avg, true_data$smoothRt_day_avg, keep = TRUE)
  
  perfect_align <- data.frame(index1 = alignment_perfect$index1, index2 = alignment_perfect$index2)
  seq_estim2 <- data.frame(index1 = alignment3$index1, index2 = alignment3$index2)
  epi_estim2 <- data.frame(index1 = alignment4$index1, index2 = alignment4$index2)
  
  dist_perfect_align <- alignment_perfect$normalizedDistance
  dist_seq_estim2 <- alignment3$normalizedDistance
  dist_epi_estim2 <- alignment4$normalizedDistance
  
  # Track which is greater and calculate proportions
  seq_greater_than_epi <- ifelse(dist_seq_estim2 > dist_epi_estim2, 1, 0)
  proportion_seq_greater <- mean(seq_greater_than_epi)
  
  # Save results for this group
  results_list[[length(results_list) + 1]] <- data.frame(
    scenario_number = scenario_number,
    dist_seq_estim2 = dist_seq_estim2,
    dist_epi_estim2 = dist_epi_estim2,
    perfect_aling= dist_perfect_align,
    proportion_seq_greater = proportion_seq_greater
  )
  
  
  # Plot histogram for outbreak data
  plot_outbreak <- ggplot(outbreak_data, aes(x = infectious_time, y = incidence)) +
    geom_bar(stat = "identity", fill = "dodgerblue4", alpha = 0.8) +
    labs(x = "Days", y = "Case Count") +
    theme_minimal()
  
  # Plot comparison of epiEstim, seqEstim, and trueRt
  plot_comparison <- ggplot() +
    geom_ribbon(data = epi_data, aes(x = t_end, ymin = `Quantile.0.05.R.`, ymax = `Quantile.0.95.R.`, fill = "casecount-Rt"), alpha = 0.3) +
    geom_line(data = epi_data, aes(x = t_end, y = `Median.R.`, color = "casecount-Rt"), size = 1) +
    geom_ribbon(data = seq_data, aes(x = time, ymin = lower, ymax = upper, fill = "sequence-Rt"), alpha = 0.4) +
    geom_line(data = seq_data, aes(x = time, y = med, color = "sequence-Rt"), size = 1) +
    geom_point(data = true_data, aes(x = infectious_time, y = smoothRt_day_avg, colour = "true-Rt"), size = 1) +
    labs(x = "Time ", y = "Rt Value") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(name = " ",
                       values = c("casecount-Rt" = "dodgerblue", "sequence-Rt" = "firebrick", "true-Rt" = "forestgreen")) +
    scale_fill_manual(name = " ", values = c("casecount-Rt" = "dodgerblue", "sequence-Rt" = "firebrick", "true-Rt" = "forestgreen"), guide = 'none')
  
  # Plot warping function alignments
  plot_warp <- ggplot() +
    geom_line(data = perfect_align, aes(x = index1, y = index2, color = "Perfect Align"), size = 1) +
    geom_line(data = seq_estim2, aes(x = index1, y = index2, color = "Sequence Estim"), size = 1) +
    geom_line(data = epi_estim2, aes(x = index1, y = index2, color = "Epi Estim"), size = 1) +
    labs(x = "Index 1 ", y = "Index 2") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(name = " ", 
                       values = c("Perfect Align" = "forestgreen", "Sequence Estim" = "firebrick", "Epi Estim" = "dodgerblue")) +
    theme(legend.position = "bottom")
  
  # Combine plots into a single row
  row_plot <- plot_outbreak + plot_comparison + plot_warp + plot_layout(ncol = 3)
  
  # Align data by time points using interpolation
  aligned_data <- data.frame(
    time = true_data$infectious_time,
    true_Rt = true_data$smoothRt_day_avg,
    seq_Rt = approx(seq_data$time, seq_data$med, xout = true_data$infectious_time, rule = 2)$y,
    epi_Rt = approx(epi_data$t_end, epi_data$Median.R., xout = true_data$infectious_time, rule = 2)$y
  )
  
  # Calculate absolute errors
  aligned_data$seq_error <- abs(aligned_data$seq_Rt - aligned_data$true_Rt)
  aligned_data$epi_error <- abs(aligned_data$epi_Rt - aligned_data$true_Rt)
  
  # Compute RMSE
  rmse_sequenceRt <- sqrt(mean(aligned_data$seq_error^2, na.rm = TRUE))
  rmse_casecountRt <- sqrt(mean(aligned_data$epi_error^2, na.rm = TRUE))
  
  # Store errors and RMSE in the results list
  results_list[[length(results_list) + 1]] <- data.frame(
    scenario_number = scenario_number,
    rmse_sequenceRt = rmse_sequenceRt,
    rmse_casecountRt = rmse_casecountRt
  )
  
  # Display results
  print("Aligned Data and Errors:")
  print(aligned_data)
  
  return(list(aligned_data = aligned_data, row_plot = row_plot, results_list = results_list))
  
}

process_all_groups <- function(groups, output_csv) {
  results_list <- list()
  
  for (group in groups) {
    result <- process_group(group, results_list)
    
    if (!is.null(result$results_list)) {
      results_list <- c(results_list, result$results_list)
    }
  }
  
  expected_columns <- c("scenario_number", "dist_seq_estim2", "dist_epi_estim2", 
                        "perfect_align", "proportion_seq_greater", "rmse_sequenceRt", "rmse_casecountRt")
  
  # Remove empty data frames before binding
  results_list <- results_list[sapply(results_list, function(x) !is.null(x) && nrow(x) > 0)]
  
  # Ensure all data frames have the same structure
  results_list <- lapply(results_list, function(df) {
    missing_cols <- setdiff(expected_columns, colnames(df))
    for (col in missing_cols) {
      df[[col]] <- NA  # Add missing columns with NA values
    }
    return(df[expected_columns])  # Reorder columns to match the expected structure
  })
  
  # Safely combine results
  results_df <- bind_rows(results_list)
  
  # Save results
  write.csv(results_df, file = output_csv, row.names = FALSE)
  
  return(results_df)
}

# Main function to process all CSV files in specified directories and save the plots
# Main function to process a limited number of CSV files
process_csv_files <- function(outbreak_dir, trueRt_dir, seqEstim_dir, epiEstim_dir, output_dir, num_groups = NULL) {
  # Create the output directory if it does not exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Read and match files
  matched_groups <- read_and_match_files(outbreak_dir, trueRt_dir, seqEstim_dir, epiEstim_dir)
  
  # If num_groups is specified, subset matched_groups
  if (!is.null(num_groups)) {
    matched_groups <- matched_groups[1:min(num_groups, length(matched_groups))]
  }
  
  # Process each group of matched files and save results to CSV
  output_csv <- paste0(output_dir, "/results.csv")
  results_df <- process_all_groups(matched_groups, output_csv)
  
  # Determine last group index 
  num_groups_to_process <- length(matched_groups)
  
  # Build plots, suppress legends for all but the last
  all_plots <- lapply(seq_along(matched_groups), function(i) {
    result <- process_group(matched_groups[[i]], list())
    if (i == num_groups_to_process) {
      # Keep legend on the last row; keep x label on last row only
      result$row_plot + labs(x = "Time")
    } else {
      # Suppress legend + remove x label on earlier rows
      (result$row_plot & theme(legend.position = "none")) +
        theme(axis.title.x = element_blank())
    }
  })
  
  
  # Arrange all plots into a single grid with 3 plots per row (one row per scenario)
  combined_plots <- wrap_plots(all_plots, ncol = 1)
  
  output_file_name <- paste0(output_dir, "/70perc_combined_plots.pdf")#update
  # Save the combined plot with optimized dimensions
  num_rows <- ceiling(length(all_plots) / 3)
  ggsave(output_file_name, plot = combined_plots, width = 15, height = 10 * num_rows, limitsize = FALSE)
}


# specify directory here
outbreak_dir <- "output/gen_incid_seventyperc_outbreaks" 
trueRt_dir <- "output/Rt_practical"
seqEstim_dir <- "output/seqRt_70perc"
epiEstim_dir <- "output/seventyperc_EPiEstim_Rt"
output_dir <- "output/figures_70perc"
process_csv_files(outbreak_dir, trueRt_dir, seqEstim_dir, epiEstim_dir, output_dir, num_groups =20)

result <- read.csv("output/figures_70perc/results.csv")


# Reshape data to long format for ggplot
result_long <- result %>%
  select(scenario_number, rmse_sequenceRt, rmse_casecountRt) %>%
  tidyr::pivot_longer(cols = c(rmse_sequenceRt, rmse_casecountRt), 
                      names_to = "RMSE_Type", values_to = "RMSE_Value")

# Define color-blind-friendly palette (Okabe-Ito: Orange & Teal)
cb_palette <- c("rmse_sequenceRt" = "#E69F00",  # Orange
                "rmse_casecountRt" = "#56B4E9") # Teal

# Create the grouped bar plot
plot_rmse <- ggplot(result_long, aes(x = scenario_number, y = RMSE_Value, fill = RMSE_Type)) +
  geom_bar(stat = "identity", position = "dodge") +  # Dodge for side-by-side bars
  labs(x = "Outbreak", y = "RMSE Value", title = " ") +
  scale_fill_manual(name = " ", values = cb_palette) +
  theme_minimal() +
  theme(legend.position = "bottom")  # Moves legend to top for better readability


# Define the file path
# file_path <- "output/figures_70perc/rmse_70perc.rds"  

#Save the ggplot object as an RDS file (without width & height)
# saveRDS(plot_rmse, file = file_path)



result_optim <- read.csv("figures_optimal/results.csv")
result_70perc <- read.csv("figures_70perc/results.csv")
result_10perc <- read.csv("figures_10perc/results.csv")
result_10to70 <- read.csv("figures_10_to_70/results.csv")
result_70to10 <- read.csv("figures_70_to_10/results.csv")

# Load data
result_optim   <- read.csv("figures_optimal/results.csv")
result_70perc  <- read.csv("figures_70perc/results.csv")
result_10perc  <- read.csv("figures_10perc/results.csv")
result_10to70  <- read.csv("figures_10_to_70/results.csv")
result_70to10  <- read.csv("figures_70_to_10/results.csv")

# Helper function to extract medians
get_medians <- function(df) {
  c(
    rmse_sequenceRt  = median(df$rmse_sequenceRt, na.rm = TRUE),
    rmse_casecountRt = median(df$rmse_casecountRt, na.rm = TRUE)
  )
}

# Create tidy data frame
scenario_names <- c("Scenario A", "Scenario B", "Scenario C", "Scenario D", "Scenario E")
all_results <- rbind(
  get_medians(result_optim),
  get_medians(result_70perc),
  get_medians(result_10perc),
  get_medians(result_10to70),
  get_medians(result_70to10)
)


df_plot <- as.data.frame(all_results) %>%
  rownames_to_column(var = "Scenario") %>%
  mutate(Scenario = scenario_names) %>%
  pivot_longer(
    cols = starts_with("rmse"),
    names_to = "Metric",
    values_to = "RMSE"
  )


cb_palette <- c(
  "rmse_sequenceRt"  = "#E69F00",  # Orange
  "rmse_casecountRt" = "#56B4E9"   # Blue
)

plot_med_rmse <- ggplot(df_plot, aes(x = Scenario, y = RMSE, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = cb_palette, name = "") +
  labs(
    title = "",
    x = "",
    y = "Median RMSE"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

# Save the ggplot object

# saveRDS(plot_med_rmse, file = "plot_med_rmse.rds")


#combine rsme plots 


# You can use `readRDS()` to load the RSD plots that were saved as R objects
rsd_optimal <- readRDS("~/Desktop/sfu-Downloads/Rt-project/elisha_siavash/figures_optimal/rmse_optim.rds")
rds_70perc <- readRDS("~/Desktop/sfu-Downloads/Rt-project/elisha_siavash/figures_70perc/rmse_70perc.rds" )
rsd_10perc <- readRDS("~/Desktop/sfu-Downloads/Rt-project/elisha_siavash/figures_10perc/rmse_10perc.rds")
rsd_10_to_70 <- readRDS("~/Desktop/sfu-Downloads/Rt-project/elisha_siavash/figures_10_to_70/rmse_10_to_70.rds")
rsd_70_to_10 <- readRDS("~/Desktop/sfu-Downloads/Rt-project/elisha_siavash/figures_70_to_10/rmse_70_to_10.rds")
rds_plot_med_rsme <- readRDS("plot_med_rmse.rds")
# Combine the plots into a grid and add labels
common_ylim <- c(0, 0.5) 

rsd_optimal   <- rsd_optimal   + ggtitle("Scenario A") + ylim(common_ylim)
rds_70perc    <- rds_70perc    + ggtitle("Scenario B") + ylim(common_ylim)
rsd_10perc    <- rsd_10perc    + ggtitle("Scenario C") + ylim(common_ylim)
rsd_10_to_70  <- rsd_10_to_70  + ggtitle("Scenario D") + ylim(common_ylim)
rsd_70_to_10  <- rsd_70_to_10  + ggtitle("Scenario E") + ylim(common_ylim)
rds_plot_med_rsme <- rds_plot_med_rsme  + ggtitle("Median RMSE") + ylim(common_ylim)

#pdf_images <- image_read_pdf("RMSE_Comparison.pdf", density = 110)  # you can adjust DPI for quality





# Combine with ggplots
combined_plot <- grid.arrange(
  rsd_optimal, rds_70perc,
  rsd_10perc, rsd_10_to_70,
  rsd_70_to_10, rds_plot_med_rsme, 
  ncol = 2, nrow = 3
)



# Save the combined plot as a PDF
ggsave("rmse_combined_plots.pdf", plot = combined_plot, width = 12, height = 10)










