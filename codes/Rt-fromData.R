#setwd("~/Downloads/Rt-project/elisha_siavash/data/data") #
library(zoo)

Inf_process_host_data_Rt_modified <- function(file_path) {
  host_data <- read.csv(file_path)
  
  # For each infector, count the number of people they infected over the entire period
  infector_counts <- host_data %>%
    group_by(infector_id) %>%
    summarise(total_infections_by_infector = n_distinct(host_id, na.rm = TRUE))
  
  # Function to calculate Rt for a single day
  calculate_day_Rt <- function(day) {
    hosts_today <- host_data %>% filter(infectious_time == day) %>% pull(host_id)
    today_infections_data <- left_join(tibble(host_id = hosts_today), infector_counts, by = c("host_id" = "infector_id"))
    
    total_infections_by_today_hosts <- sum(today_infections_data$total_infections_by_infector, na.rm = TRUE)
    avg_infections_per_host_today <- total_infections_by_today_hosts / length(hosts_today)
    
    return(avg_infections_per_host_today)
  }
  
  # Apply the function for each day
  Rt_fromData <- host_data %>%
    distinct(infectious_time) %>%
    rowwise() %>%
    mutate(Rt_day_avg = calculate_day_Rt(infectious_time))
  
  # Ensure all days are represented
  all_days <- tibble(infectious_time = seq(min(host_data$infectious_time), max(host_data$infectious_time), by=1))
  Rt_fromData <- left_join(all_days, Rt_fromData, by="infectious_time")
  
  # Fill NAs for days without exposures
  Rt_fromData$Rt_day_avg[is.na(Rt_fromData$Rt_day_avg)] <- 0
  
  return(Rt_fromData)
}






Inf_process_host_data_Rt_modified_weekly <- function(file_path) {
  host_data <- read.csv(file_path)
  
  # Convert infectious_time to a weekly period
  host_data <- host_data %>%
    mutate(week_period = ceiling(infectious_time / 7))
  
  # For each infector, count the number of people they infected over the entire period
  infector_counts <- host_data %>%
    group_by(infector_id) %>%
    summarise(total_infections_by_infector = n_distinct(host_id, na.rm = TRUE))
  
  # Function to calculate Rt for a single week
  calculate_week_Rt <- function(week) {
    hosts_this_week <- host_data %>% 
      filter(week_period == week) %>% 
      pull(host_id)
    this_week_infections_data <- left_join(tibble(host_id = hosts_this_week), infector_counts, by = c("host_id" = "infector_id"))
    
    total_infections_by_this_week_hosts <- sum(this_week_infections_data$total_infections_by_infector, na.rm = TRUE)
    avg_infections_per_host_this_week <- total_infections_by_this_week_hosts / length(hosts_this_week)
    
    return(avg_infections_per_host_this_week)
  }
  
  # Apply the function for each week
  Rt_fromData_weekly <- host_data %>%
    distinct(week_period) %>%
    rowwise() %>%
    mutate(Rt_week_avg = calculate_week_Rt(week_period))
  
  # Ensure all weeks are represented
  all_weeks <- tibble(week_period = seq(min(host_data$week_period), max(host_data$week_period), by=1))
  Rt_fromData_weekly <- left_join(all_weeks, Rt_fromData_weekly, by="week_period")
  
  # Fill NAs for weeks without exposures
  Rt_fromData_weekly$Rt_week_avg[is.na(Rt_fromData_weekly$Rt_week_avg)] <- 0
  
  return(Rt_fromData_weekly)
}
