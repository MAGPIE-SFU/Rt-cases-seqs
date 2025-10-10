library(dplyr)
library(EpiEstim)

# Function to handle estimation with retries
estimate_R_with_retry <- function(incid_data, gen_time_data, max_retries = 10) {
  for (retry in 0:max_retries) {
    t_start <- seq(2, nrow(incid_data) - (7 + retry))
    t_end <- t_start + (7 + retry)
    
    result <- tryCatch({
      Rt_res <- estimate_R(incid = incid_data$incidence, 
                           method = "parametric_si",
                           config = make_config(list(
                             mean_si = mean(gen_time_data$generation_times), 
                             std_si = sd(gen_time_data$generation_times), 
                             t_start = t_start, 
                             t_end = t_end)))
      list(result = Rt_res, warning = NULL)
    }, warning = function(w) {
      message("Warning in estimate_R: ", w$message, " - Retrying with adjusted t_start and t_end...")
      list(result = NULL, warning = w)
    }, error = function(e) {
      message("Error in estimate_R: ", e$message)
      stop("Failed due to error.")
    })
    
    if (is.null(result$warning)) {
      return(result$result)
    }
  }
  stop("All retries failed")
}

# Set the input and output directories
input_dir <- "output/gen_incid_tenperc_outbreaks"
# set your own output folder. If you uncomment the next line, you will re-write this one
# output_dir <- "output/tenperc_EPiEstim_Rt"

# Create the output directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Get a list of all .csv files in the input directory
gen_files <- list.files(input_dir, pattern = "*_gt.csv", full.names = TRUE)
incid_files <- list.files(input_dir, pattern = "*_incid.csv", full.names = TRUE)

# Loop over each file, apply the function, and save the result

# NOTE-- this will write in the output folder! 
for (i in seq_along(incid_files)) {
  incid_dat_Rt <- read.csv(incid_files[i])
  gen_time <- read.csv(gen_files[i])
  
  # Try to estimate R and catch errors to skip over problematic files
  tryCatch({
    Rt_res <- estimate_R_with_retry(incid_dat_Rt, gen_time)
    
    # Save the Rt results
    base_name <- tools::file_path_sans_ext(basename(incid_files[i]))
    output_file_name <- file.path(output_dir, paste0(base_name, "_EpiEstimRt.csv"))
    
    Rt_data <- as.data.frame(Rt_res$R)
    write.csv(Rt_data, file = output_file_name, row.names = FALSE)
  }, error = function(e) {
    message("Skipping file due to repeated errors: ", incid_files[i])
  })
}
