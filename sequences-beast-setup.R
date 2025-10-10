# Loading the packges
library(Biostrings)
library(readr)
library(dplyr)
library(xml2)

# This script brings in the sequence data from a simulated outbreak (in fasta format)
# together with additional information (e.g. name, host, collection time). 
# It then has several functions to sample the outbreak and return a downsampled
# fasta file and corresponding host information, along with an xml, with which
# Rt can be estimated in BEAST. 


# Reading the files
directory <- "data_outbreaks-diff"
fasta_file_name <- "outbreak_001_wgs"
meta_data <- "outbreak_001_wgs_metadata"
xml_template_file <- file.path("template.xml")
setwd(directory)

# Read the FASTA file
fasta_file <- readDNAStringSet(paste(fasta_file_name,".fasta",sep = ""))

# Read the CSV file with sample IDs
meta_data_df <- read_csv(paste(meta_data,".csv",sep=""))

# Create a dataframe from meta data (all data) that shows the cumulative recurrence in each time 
meta_check <- meta_data_df %>%
  group_by(collection_time) %>%
  summarize(recurrence = n())

# calculating proportion in a way that sum of sequences become 400
max_seq <- 400

# basic proportion of sampling
pro1 <- round(max_seq/sum(meta_check$recurrence),2)

# subsampling rate 1
pro2 <- pro1*0.7

# subsampling rate 2
pro3 <- pro1*0.1

# calculating limit for the sampling with limit
lim <- round(1/pro1)

# finding time of the peak
peak_time <- min(meta_check$collection_time[meta_check$recurrence==max(meta_check$recurrence)])
plot(recurrence~collection_time,meta_check,type='l',xlab="time")
abline(v = peak_time, col = "red", lwd = 2) 

# Below are three different functions for sampling. 
# Function 1: 1 proportion sampling
# Function 2: Optimum case, 
# Function 3: 2 Proportions sampling

# Function 1: One proportion sampling regime
prop_sampling <- function(prop){
  dir_name <- paste("prop_",prop,sep = "")
  if (!dir.exists(dir_name)) {
    dir.create(dir_name, recursive = TRUE)
    
  }
  setwd(dir_name)
  # creating the meta data of sampled population 
  sampled_meta <- meta_data_df %>% sample_frac(prop)
  meta_file_name1 <- paste(meta_data,'_prop_',prop,sep = "")
  meta_file_name2 <- paste(meta_file_name1,".csv",sep = "")
  
  # writing the meta data of sampled recurrence in a file 
  write.csv(sampled_meta,meta_file_name2)
  
  # Create the new dataframe to calculate the summary of recurrence in sampled 
  sampled_check <- sampled_meta %>%
    group_by(collection_time) %>%
    summarize(recurrence = n())
  
  # writing the summary of sampled recurrence in a file 
  write.csv(sampled_check,paste(meta_data,'_samp_recur_prop_',prop,".csv",sep = ""))
  
  plot_title <- paste("prop ",prop,",total ",sum(sampled_check$recurrence))
  png(paste(plot_title,".png",sep = ""), width = 800, height = 600, res = 100)
  plot(recurrence~collection_time,meta_check,type='l',xlab="time")
  points(recurrence~collection_time,sampled_check)
  title(main = plot_title)
  dev.off()
  
  # Extract the IDs from the 'name' column
  ids_to_extract <- sampled_meta$name # in Rt project run this line
  sampled_fasta <- fasta_file[names(fasta_file) %in% ids_to_extract] 
  
  fasta_samp_file_name1 <- paste(fasta_file_name,"_prop_",prop,sep = "")
  fasta_samp_file_name2 <- paste(fasta_samp_file_name1,".fasta",sep = "")
  
  # Write the filtered sequences to a new FASTA file
  writeXStringSet(sampled_fasta, fasta_samp_file_name2)
  
  return(c(meta_file_name1,fasta_samp_file_name1))
  
}


# Sampling based on rate (basic) and a limit (optimum case)
prop_sampling_lim <- function(prop,lim){
  dir_name <- paste("opt_prop_",prop,"_lim_",lim,sep = "")
  if (!dir.exists(dir_name)) {
    dir.create(dir_name, recursive = TRUE)
    
  }
  setwd(dir_name)
  # Join recurrence back to the original dataframe
  meta_data_with_recurrence <- meta_data_df %>%
    left_join(meta_check, by = "collection_time")
  
  # Sample based on the criteria
  sampled_meta <- meta_data_with_recurrence %>%
    group_by(collection_time) %>%
    mutate(sample_indicator = case_when(
      recurrence > lim ~ row_number() %in% sample(row_number(), size = max(1, round(recurrence * pro1))),
      TRUE ~ row_number() == sample(row_number(), size = 1)
    )) %>%
    filter(sample_indicator) %>%
    ungroup() %>%
    select(-recurrence, -sample_indicator)
  
  
  meta_file_name1 <- paste(meta_data,"_",dir_name,sep = "")
  meta_file_name2 <- paste(meta_file_name1,".csv",sep = "")
  
  # writing the meta data of sampled recurrence in a file 
  write.csv(sampled_meta,meta_file_name2)
  
  # Create the new dataframe to calculate the summary of recurrence in sampled 
  sampled_check <- sampled_meta %>%
    group_by(collection_time) %>%
    summarize(recurrence = n())
  
  # writing the summary of sampled recurrence in a file 
  write.csv(sampled_check,paste(meta_data,'_samp_recur_',dir_name,".csv",sep = ""))
  
  plot_title <- paste("prop ",prop,',lim ',lim,",total ",sum(sampled_check$recurrence),sep = "")
  png(paste(plot_title,".png",sep = ""), width = 800, height = 600, res = 100)
  plot(recurrence~collection_time,meta_check,type='l',xlab="time")
  points(recurrence~collection_time,sampled_check)
  title(main = plot_title)
  dev.off()
  
  
  # Extract the IDs from the 'name' column
  ids_to_extract <- sampled_meta$name # in Rt project run this line
  sampled_fasta <- fasta_file[names(fasta_file) %in% ids_to_extract] 
  
  fasta_samp_file_name1 <- paste(fasta_file_name,"_",dir_name,sep = "")
  fasta_samp_file_name2 <- paste(fasta_samp_file_name1,".fasta",sep = "")
  
  # Write the filtered sequences to a new FASTA file
  writeXStringSet(sampled_fasta, fasta_samp_file_name2)
  
  return(c(meta_file_name1,fasta_samp_file_name1))
  
}


# Function 3: Two proportions sampling
two_prop <- function(prop2,prop3){
  # Create sampled_meta3 with two sampling proportions
  # the idea is to change the proportion a bit before the peak time
  dir_name <- paste("prop_",prop2,"_prop_",prop3,sep = "")
  if (!dir.exists(dir_name)) {
    dir.create(dir_name, recursive = TRUE)
    
  }
  setwd(dir_name)
  
  peak_time <- peak_time-10
  sampled_meta3 <- meta_data_df %>%
    mutate(
      sampled = if_else(collection_time <= peak_time, row_number() %in% sample(row_number(), n() * prop2), row_number() %in% sample(row_number(), n() * prop3))
    ) %>%
    filter(sampled) %>%
    select(-sampled)
  
  # Create the new dataframe to calculate the summary of recurrence in sampled 
  sampled_check3 <- sampled_meta3 %>%
    group_by(collection_time) %>%
    summarize(recurrence = n())
  
  plot_title <-  paste("cut-off time ",peak_time,", ","rates ",prop2,", ",prop3,"total ",sum(sampled_check3$recurrence),sep = "")
  png(paste(plot_title,".png",sep = ""), width = 800, height = 600, res = 100)
  plot(recurrence~collection_time,meta_check,type='l',xlab="time")
  points(recurrence~collection_time,sampled_check3)
  title(main =plot_title)
  dev.off()
  
  # Extract the IDs from the 'name' column
  ids_to_extract3 <- sampled_meta3$name
  sampled_fasta3 <- fasta_file[names(fasta_file) %in% ids_to_extract3] 
  
  meta_file_name1 <- paste(meta_data,'_time_',peak_time,"_rates_",prop2,"_",prop3,sep = "")
  meta_file_name2 <- paste(meta_file_name1,".csv",sep = "")
  
  # writing the meta data of sampled recurrence in a file 
  write.csv(sampled_meta3,meta_file_name2)
  
  # writing the summary of sampled recurrence in a file 
  write.csv(sampled_check3,paste(meta_data,"_samp_recur_time_",peak_time,"_rates_",prop2,"_",prop3,".csv",sep = ""))
  
  fasta_samp_file_name1 <- paste(fasta_file_name,'_time_',peak_time,"_prop_",prop2,"_",prop3,sep = "")
  fasta_samp_file_name2 <- paste(fasta_samp_file_name1,".fasta",sep = "")
  
  # Write the filtered sequences to a new FASTA file
  writeXStringSet(sampled_fasta3, fasta_samp_file_name2)
  
  return(c(meta_file_name1,fasta_samp_file_name1))
}

# Adding date to fasta file from meta data
add_date <- function(raw_file_names){
  meta_file_name <- raw_file_names[1]
  fasta_file_name <- raw_file_names[2]
  output_fasta_name1 = paste(fasta_file_name,"_full",sep = "")
  output_fasta_name2 = paste(output_fasta_name1,".fasta",sep = "")
  
  meta = read.csv(paste(meta_file_name,".csv",sep = ""))

  name_to_time <- setNames(meta$collection_time, meta$name)
  
  # Read the fasta file
  fasta_seq <- readDNAStringSet(paste(fasta_file_name,".fasta",sep = ""))
  
  # Modify the fasta headers
  names(fasta_seq) <- sapply(names(fasta_seq), function(name) {
    paste(name, name_to_time[name], sep="|")
  })
  
  # Write the modified fasta to a new file
  writeXStringSet(fasta_seq, output_fasta_name2)
  return(c(output_fasta_name1,output_fasta_name2))
  
}

# Chaning xml template file according to updated fasta file
# Function to read a file as a single string
read_file_as_string <- function(file_path) {
  lines <- readLines(file_path)
  fasta_string <- paste(lines, collapse = "\n")
  return(fasta_string)
}

# Function to parse fasta string and return sequence IDs, taxon values, dates, and sequence values
parse_fasta <- function(fasta_string) {
  # Splitting the FASTA string into individual lines
  lines <- unlist(strsplit(fasta_string, "\n"))
  
  # Initialize lists to store the extracted data
  ids <- c()
  dates <- c()
  sequences <- c()
  
  # Temporary variable to accumulate sequence data
  current_sequence <- ""
  for (line in lines) {
    if (startsWith(line, ">")) {
      # Save the current sequence (if it exists) before starting a new one
      if (nchar(current_sequence) > 0) {
        sequences <- c(sequences, current_sequence)
        current_sequence <- "" # Reset sequence accumulator
      }
      # Extract ID and additional info
      id_info <- strsplit(sub("^>", "", line), "\\|")[[1]]
      ids <- c(ids, id_info[1])
      # Safely extract taxon and date if they exist
      if (length(id_info) >= 2) {
        dates <- c(dates, id_info[2])
      } else {
        dates <- c(dates, NA) # Use NA for missing data
      }
    } else {
      current_sequence <- paste0(current_sequence, line)
    }
  }
  
  # Don't forget to add the last sequence
  if (nchar(current_sequence) > 0) {
    sequences <- c(sequences, current_sequence)
  }
  
  return(list(ids = ids, dates = dates, values = sequences))
}

# Function to update XML with fasta data
update_xml_with_fasta <- function(fasta_file_input, xml_file_input, xml_file_ouput) {
  
  # Read fasta file as string and parse data
  fasta_string <- read_file_as_string(fasta_file_input)
  fasta_data <- parse_fasta(fasta_string)
  # Load the XML file
  doc <- read_xml(xml_file_input)
  
  # Update the <data> id attribute to 'file_name'
  data_node <- xml_find_first(doc, "//data")
  
  # Remove all existing sequence nodes
  xml_find_all(doc, "//sequence") %>% xml_remove()
  
  # Add new sequence nodes based on fasta data
  for (i in seq_along(fasta_data$ids)) {
    xml_add_child(data_node, "sequence", 
                  id = paste0(fasta_data$ids[i], "|", fasta_data$dates[i]), 
                  spec = "Sequence", 
                  taxon = fasta_data$ids[i], 
                  totalcount = "4", 
                  value = fasta_data$values[i])
  }
  
  # Construct the value string for the trait node
  trait_value <- paste(paste(fasta_data$ids, fasta_data$dates, sep="="), collapse=",")
  
  # Find the trait node and update its value
  trait_node <- xml_find_first(doc, "//trait[@id='dateTrait.t:template']")
  xml_set_attr(trait_node, "value", trait_value)
  
  # Save the updated XML
  write_xml(doc,paste(xml_file_ouput,".xml",sep = "") )  # Replace with your desired output file path
}

# High proportion of sampling
file_prop <- prop_sampling(pro2)
dated_prop <- add_date(file_prop)
update_xml_with_fasta(dated_prop[2], xml_template_file, file_prop[2])
setwd(directory)

# Low proportion of sampling
file_prop <- prop_sampling(pro3)
dated_prop <- add_date(file_prop)
update_xml_with_fasta(dated_prop[2], xml_template_file, file_prop[2])
setwd(directory)

# One proportion and a limit of sampling (Optimum)
file_prop_lim <- prop_sampling_lim(pro1,lim)
dated_prop_lim <- add_date(file_prop_lim)
update_xml_with_fasta(dated_prop_lim[2], xml_template_file, file_prop_lim[2])
setwd(directory)

# From high proportion to low
file_high_low <- two_prop(pro2,pro3)
dated_high_low <- add_date(file_high_low)
update_xml_with_fasta(dated_high_low[2], xml_template_file, file_high_low[2])
setwd(directory)

# From low proportion to high
file_low_high <- two_prop(pro3,pro2)
dated_low_high <- add_date(file_low_high)
update_xml_with_fasta(dated_low_high[2], xml_template_file, file_low_high[2])
setwd(directory)

