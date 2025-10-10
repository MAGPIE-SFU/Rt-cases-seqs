
install.packages("remote")
library(remote)

# install.packages("devtools")
library(devtools)

install_github("jessicastockdale/OOPidemic")

library(OOPidemic)

set.seed(543)

if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
}



if (!dir.exists("data_outbreaks")) {
  dir.create("data_outbreaks")
}


num_outbreaks <- 20 # Number of outbreaks

# NOTE that this loop writes to the folder data_outbreaks-diff
for (i in 1:num_outbreaks) {
  # Sample total population and initial infected from uniform distribution
  total_pop <- round(runif(1, min = 10000, max = 15000)) # Between 10000 and 15000
  init_inf <- round(runif(1, min = 5, max = 10)) # Between 5 and 10
  inf_rate <- runif(1, min = 0.21, max = 0.25) # Between 0.1 and 0.25
  # Initialize the reference strain
  ref_strain <- ReferenceStrain$new(
    name = "ref_strain",
    g_len = 1000,
    mut_rate = 1 / 700 # Mutation rate
  )
  
  init_sus <- total_pop - init_inf
  
  # Initialize the group with the new total_pop, init_inf, and varying inf_rate
  group <- Group$new(
    id = 1,
    ref_strain = ref_strain,
    init_inf = init_inf,
    init_sus = init_sus,
    max_init_dist = 3,
    inf_rate = inf_rate, # Varying infection rate
    rec_shape = 20, rec_rate = 4 # Recovery settings,  
    #rec_shape-=25 increases duration
  )
  
  lab <- Lab$new()
  
  # Simulate the outbreak
  index_cases <- group$infectious_hosts()
  any_case <- group$hosts
  
  lab$sample_hosts(group$hosts_due_for_sampling, group$time)
  
  while (all(
    group$is_outbreak_active,
    group$recovered_size < total_pop, # Adjust according to total_pop
    group$time < 100000
  )) {
    group$infect()
    lab$sample_hosts(group$hosts_due_for_sampling, group$time)
  }
  print(c(total_pop,init_inf,inf_rate ))
  # Generate the filename for the current outbreak
  filename <- sprintf("outbreak_%03d", i) # Fixed spacing in format string
  
  # Save the outbreak data in the 'data_outbreaks' folder
  lab$save("data_outbreaks-diff", filename) #older version data_outbreaks
}


