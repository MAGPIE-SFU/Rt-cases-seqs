# Load required libraries
library(MASS)
library(dplyr)

# this file has a downsampling function and a function to process input csv files.
# Be careful with the output folders as it would be relatively easy to over-write 
# the current output. Relevant commands are commented out (we think). 

# downsamples the line list
sample_line_list <- function(data, size, prob) {
  # Function to perform the sampling
  sample_ids <- function(data, size, prob) {
    group_size <- nrow(data)
    sample_size <- rnbinom(1, size = size, prob = prob)
    sample_size <- min(sample_size, group_size) # Ensuring sample size is not larger than the group
    sampled_data <- data[sample(1:group_size, size = sample_size), ]
    return(sampled_data)
  }
  
  # Group by exposure_time and apply sampling
  sampled_line_list <- data %>%
    group_by(exposure_time) %>%
    do(sample_ids(., size, prob)) %>%
    ungroup() # Make sure to ungroup to avoid grouping issues in later steps
  
  return(sampled_line_list)
}

# Main function to process all CSV files in the input directory
process_csv_files <- function(input_dir, output_dir, fraction_to_keep) {
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
      mean_downsampled <- group_sizes$n * fraction_to_keep
      size_estimate <- mean(group_sizes$n, na.rm = TRUE)
      prob_downsampled <- size_estimate / (size_estimate + mean(mean_downsampled, na.rm = TRUE))
      
      # Check for NA values
      if (is.na(size_estimate) || is.na(prob_downsampled)) {
        warning(paste("NA values found in parameter estimation for file:", file_path))
        next
      }
      
      # Downsample the line list
      sampled_data <- sample_line_list(data, size_estimate, prob_downsampled)
      
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

# Specify dir
# for example: 
input_dir <- "data_outbreaks-diff"
output_dir <- "output/tenperc_outbreaks_downsampled"
fraction_to_keep <- 0.1

# example of usage: 
# process_csv_files(input_dir, output_dir, fraction_to_keep)
