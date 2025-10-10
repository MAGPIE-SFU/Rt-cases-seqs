# Load necessary library
library(dplyr)



# Function to modify a dataframe, creating a rolling average 
modify_dataframe <- function(seqRt_data) {
  colnames(seqRt_data) <- NULL
  
  transp_seqRt_data <- as.data.frame(t(seqRt_data))
  
  colnames(transp_seqRt_data) <- unlist(transp_seqRt_data[1, ])
  transp_seqRt_data <- transp_seqRt_data[-1, ]
  transp_seqRt_data <- transp_seqRt_data %>% mutate(time = 1:nrow(transp_seqRt_data))
  
  transp_seqRt_data[] <- lapply(transp_seqRt_data, as.numeric)
  
  transp_seqRt_data <- transp_seqRt_data[rev(rownames(transp_seqRt_data)), ]
  transp_seqRt_data$time <- 1:nrow(transp_seqRt_data)
  
  # Create a new dataframe 'ave_transp_seqRt_data' with the rolling average
  ave_transp_seqRt_data <- data.frame(
    lower = sapply(seq(1, nrow(transp_seqRt_data)-2), function(i) mean(transp_seqRt_data$lower[i:(i+2)], na.rm = TRUE)),
    med = sapply(seq(1, nrow(transp_seqRt_data)-2), function(i) mean(transp_seqRt_data$med[i:(i+2)], na.rm = TRUE)),
    upper = sapply(seq(1, nrow(transp_seqRt_data)-2), function(i) mean(transp_seqRt_data$upper[i:(i+2)], na.rm = TRUE)),
    time = sapply(seq(1, nrow(transp_seqRt_data)-2), function(i) mean(transp_seqRt_data$time[i:(i+2)], na.rm = TRUE))
  )
  
  return(ave_transp_seqRt_data)
}

# Function to process all CSV files in a directory and save the modified data
process_csv_files <- function(input_dir, output_dir) {
  # Create the output directory if it does not exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
    print(paste("Created directory:", output_dir))
  }
  
  # Get a list of all .csv files in the input directory
  file_list <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
  print(paste("Files found:", length(file_list)))
  
  # Loop over each file, apply the function, and save the result
  for (file_path in file_list) {
    print(paste("Processing file:", file_path))
    seqRt_data <- read.csv(file_path)
    
    # Check if the data is read correctly
    print(head(seqRt_data))
    
    modified_data <- modify_dataframe(seqRt_data)
    
    # Create the output file name
    base_name <- tools::file_path_sans_ext(basename(file_path))
    output_file_name <- paste0(output_dir, "/", base_name, "_modified.csv")
    
    # Save the result to a new .csv file
    write.csv(modified_data, file = output_file_name, row.names = FALSE)
    print(paste("Saved modified file to:", output_file_name))
  }
}

# Specify dir
input_dir <- "~/Downloads/Rt-project/final-analysis/70_10"
output_dir <- "~/Downloads/Rt-project/elisha_siavash/seqRt_70_to_10"
process_csv_files(input_dir, output_dir)
