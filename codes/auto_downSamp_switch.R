library(MASS)
library(dplyr)

# Function to perform the sampling within each group
sample_ids <- function(data, size, prob) {
  group_size <- nrow(data)
  sample_size <- rnbinom(1, size = size, prob = prob)
  sample_size <- min(sample_size, group_size) # Ensuring sample size is not larger than the group
  sampled_data <- data[sample(1:group_size, size = sample_size), ]
  return(sampled_data)
}

# Function to downsample a line list with switching fractions
sample_line_list <- function(data, size, prob, switch_time, new_size, new_prob) {
  sampled_data <- data %>%
    group_by(exposure_time) %>%
    do({
      if (unique(.$exposure_time) < switch_time) {
        sample_ids(., size, prob)
      } else {
        sample_ids(., new_size, new_prob)
      }
    }) %>%
    ungroup() # Make sure to ungroup to avoid grouping issues in later steps
  
  return(sampled_data)
}

# Main function to process all CSV files in the input directory
process_csv_files <- function(input_dir, output_dir, initial_fraction_to_keep, switch_time, new_fraction_to_keep) {
  # Create the output directory if it does not exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Get a list of all .csv files in the input directory with the specified naming convention
  file_list <- list.files(input_dir, pattern = "_hostdata.*\\.csv$", full.names = TRUE)
  
  # Process each file
  for (file_path in file_list) {
    tryCatch({
      # Read the data
      data <- read.csv(file_path)
      
      # Ensure the necessary columns are present
      if (!all(c("host_id", "group", "infector_id", "infector_group", "exposure_time", "infectious_time") %in% colnames(data))) {
        warning(paste("File skipped due to missing columns:", file_path))
        next
      }
      
      # Estimate parameters of the negative binomial distribution from the original data
      group_sizes <- data %>% group_by(exposure_time) %>% summarise(n = n())
      mean_downsampled_initial <- group_sizes$n * initial_fraction_to_keep
      size_estimate_initial <- mean(group_sizes$n, na.rm = TRUE)
      prob_downsampled_initial <- size_estimate_initial / (size_estimate_initial + mean(mean_downsampled_initial, na.rm = TRUE))
      
      mean_downsampled_new <- group_sizes$n * new_fraction_to_keep
      size_estimate_new <- mean(group_sizes$n, na.rm = TRUE)
      prob_downsampled_new <- size_estimate_new / (size_estimate_new + mean(mean_downsampled_new, na.rm = TRUE))
      
      # Check for NA values
      if (is.na(size_estimate_initial) || is.na(prob_downsampled_initial) || is.na(size_estimate_new) || is.na(prob_downsampled_new)) {
        warning(paste("NA values found in parameter estimation for file:", file_path))
        next
      }
      
      # Downsample the line list
      sampled_data <- sample_line_list(data, size_estimate_initial, prob_downsampled_initial, switch_time, size_estimate_new, prob_downsampled_new)
      
      # Create the output file name
      base_name <- tools::file_path_sans_ext(basename(file_path))
      output_file_name <- paste0(output_dir, "/", base_name, "_downsampled.csv")
      
      # Save the downsampled data to a new .csv file
      write.csv(sampled_data, file = output_file_name, row.names = FALSE)
      print(paste("Saved downsampled file to:", output_file_name))
    }, error = function(e) {
      warning(paste("Error processing file:", file_path, "Error message:", e$message))
      next
    })
  }
}

# Specify directory 
# example: 
input_dir <- "data_outbreaks-diff"
output_dir <- "output/10_to_70_outbreaks_downsampled"
initial_fraction_to_keep <- 0.1
switch_time <- 50 # Switch after this time step
new_fraction_to_keep <- 0.7
# process_csv_files(input_dir, output_dir, initial_fraction_to_keep, switch_time, new_fraction_to_keep)
