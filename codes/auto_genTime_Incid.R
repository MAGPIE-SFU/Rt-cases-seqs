source("./data/generate_incid.R")
source("./gen_interval.R")
setwd("~/Downloads/Rt-project/elisha_siavash")



input_dir <- "output/tenperc_outbreaks_downsampled"
output_dir <- "output/gen_incid_tenperc_outbreaks"

# Create the output directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Get a list of all .csv files in the input directory
file_list <- list.files(input_dir, pattern = "*hostdata_downsampled.csv", full.names = TRUE)

# Loop over each file, apply the function, and save the result
for (file_path in file_list) {
  tryCatch({
    # Apply the calculate_generation_times function
    generation_times <- calculate_generation_times(file_path)
    
    base_name <- tools::file_path_sans_ext(basename(file_path))
    output_file_name <- paste0(output_dir, "/", base_name, "_gt.csv")
    
    # Save the result to a new .csv file
    write.csv(generation_times, file = output_file_name, row.names = FALSE)
    print(paste("Saved generation times to:", output_file_name))
  }, error = function(e) {
    warning(paste("Error processing file:", file_path, "Error message:", e$message))
  })
}

# Get a list of all .csv files in the input directory
file_list_incid <- list.files(input_dir, pattern = "*hostdata_downsampled.csv", full.names = TRUE)

# Loop over each file, apply the function, and save the result
for (file_path in file_list_incid) {
  tryCatch({
    # Apply the generate_incidence_data function
    incid_dat <- generate_incidence_data(file_path)
    
    base_name <- tools::file_path_sans_ext(basename(file_path))
    output_file_name <- paste0(output_dir, "/", base_name, "_incid.csv")
    
    # Save the result to a new .csv file
    write.csv(incid_dat, file = output_file_name, row.names = FALSE)
    print(paste("Saved incidence data to:", output_file_name))
  }, error = function(e) {
    warning(paste("Error processing file:", file_path, "Error message:", e$message))
  })
}




