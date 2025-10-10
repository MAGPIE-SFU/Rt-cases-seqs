###for methods Low-pass Filter: A low-pass filter can help remove high-frequency noise.
###library(pracma)
#smoothed_rt <- movavg(rt_data, n = 7, type = "s")
library(pracma)
source("data/Rt-fromData.R")



# Function to calculate Rt and apply smoothing with dynamic window adjustment
calculate_and_smooth_Rt <- function(file_path) {
  # Apply the calculate_generation_times function
  pract_Rt <- Inf_process_host_data_Rt_modified(file_path)
  
  # Initialize smoothing window size
  n <- 7
  smoothed_Rt <- NULL
  
  # Applying the smoothing function with decreasing window size
  while (n >= 2) {
    try({
      smoothed_Rt <- movavg(pract_Rt$Rt_day_avg, n = n, type = "s")
      break  # Exit loop if smoothing is successful
    }, silent = TRUE)
    n <- n - 1
  }
  
  if (is.null(smoothed_Rt)) {
    warning(paste("Smoothing failed for file:", file_path))
    pract_Rt$smoothRt_day_avg <- NA  # Assign NA if smoothing fails
  } else {
    pract_Rt$smoothRt_day_avg <- smoothed_Rt
  }
  
  return(pract_Rt)
}

# Set the input and output directories
input_dir <- "output/thirtyperc_outbreaks_downsampled"
output_dir <- "output/thirtyperc_Rt_practical"

# Create the output directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Get a list of all .csv files in the input directory
file_list <- list.files(input_dir, pattern = "*hostdata_downsampled.csv", full.names = TRUE)

# Loop over each file, apply the function, and save the result
for (file_path in file_list) {
  # Calculate Rt and apply smoothing
  pract_Rt <- calculate_and_smooth_Rt(file_path)
  
  # Create the output file name
  base_name <- tools::file_path_sans_ext(basename(file_path))
  output_file_name <- paste0(output_dir, "/", base_name, "_practRt.csv")
  
  # Save the result to a new .csv file
  write.csv(pract_Rt, file = output_file_name, row.names = FALSE)
}




##outbreak_1_test_daily <- Inf_process_host_data_Rt_modified("~/Downloads/Rt-project/elisha_siavash/data_outbreaks/outbreak_014_hostdata.csv")

#epiRTtest <- read.csv("EPiEstim_Rt/outbreak_014_hostdata_incid_EpiEstimRt.csv")
plot(pract_Rt$Rt_day_avg)
lines(pract_Rt$smoothRt_day_avg, type ="p" , col="red")
#lines(filtered_rt,type ="p" , col="green")
#lines(epiRTtest$Median.R., type ="p" , col="blue")




