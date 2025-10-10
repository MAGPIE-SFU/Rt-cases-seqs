# Load the required packages and set global options 
library(coda)
library(bdskytools)
library(beastio)
library(RColorBrewer)

# I believe this file collects the BEAST outputs and generates csv files like outbreak_001_wgs_time_68_prop_0.021_0.049_modified.csv
# which contain Rt estimates from the simulated sequences 


maxtime <- 154 # Last sampled time, could be found in _full.fasta file 
folder <- 8  # folder name
pro1 <- 0.009 # subsampling rate 1
pro2 <- 0.063 # subsampling rate 2
peaktime <- 64 # tiem of the pick
lim <- 11 # Limit = 1/(basic proportion of sampling) (calculated in Code 1)

D3 <- function(number) {
  # Format the number to always have three digits
  formatted_number <- sprintf("%03d", number)
  return(formatted_number)
}

input_dir <- "data_outbreaks-diff/"

# in case of one proportion
folder_dir <- paste(input_dir,folder,"/prop_",pro1,"/",sep = "")
file_name <- paste("outbreak_",D3(folder),"_wgs_prop_",pro1,sep = "")

# in case of high - low
folder_dir <- paste(input_dir,folder,"/prop_",pro1,"_prop_",pro2,"/",sep = "")
file_name <- paste("outbreak_",D3(folder),"_wgs_time_",peaktime,"_prop_",pro1,"_",pro2,sep = "")


# in case of low - high
folder_dir <- paste(input_dir,folder,"/prop_",pro2,"_prop_",pro1,"/",sep = "")
file_name <- paste("outbreak_",D3(folder),"_wgs_time_",peaktime,"_prop_",pro2,"_",pro1,sep = "")



# in case of optimum (proportion + limit)
folder_dir <- paste(input_dir,folder,"/opt_prop_",pro1,"_lim_",lim,"/",sep = "")
file_name <- paste("outbreak_",D3(folder),"_wgs_opt_prop_",pro1,"_lim_",lim,sep = "")

file_dir <- paste(folder_dir,file_name,".log",sep = "")
setwd(folder_dir)

# Set up colours
cols  <- list(blue   = RColorBrewer::brewer.pal(12,"Paired")[2], 
              orange = RColorBrewer::brewer.pal(12,"Paired")[8])

set_alpha <- function(c, alpha=1.0) paste0(c,format(as.hexmode(round(alpha*255)), width=2))

# Parameters
params <- list(
  logfile = file_dir,
  gridsize = maxtime,
  mostrecent = "2020-12-29"
)

mostrecent_decimal <- lubridate::decimal_date(lubridate::ymd(params$mostrecent))

# Load the trace file
bdsky_trace   <- beastio::readLog(params$logfile, burnin=0.1)

# Check convergence
beastio::checkESS(bdsky_trace)

# Plot the ESS values of all parameters
beastio::checkESS(bdsky_trace, cutoff=200, plot=TRUE, log='y', ylim=c(1,10000), title="All parameters", plot.grid=TRUE)

# Extract parameter estimates and HPDs
Re_sky <- beastio::getLogFileSubset(bdsky_trace, "reproductiveNumber_BDSKY_Serial")
Re_hpd <- t(beastio::getHPDMedian(Re_sky))
delta_hpd <- beastio::getHPDMedian(bdsky_trace[, "becomeUninfectiousRate_BDSKY_Serial"])

# Plotting non-gridded BDSKY estimates
bdskytools::plotSkyline(1:10, Re_hpd, type='step', ylab="Re")

# Gridding and plotting smooth skyline
tmrca_med  <- median(bdsky_trace[, "origin_BDSKY_Serial"])
gridTimes  <- seq(0, median(tmrca_med), length.out=params$gridsize)  
Re_gridded <- mcmc(bdskytools::gridSkyline(Re_sky, bdsky_trace[, "origin_BDSKY_Serial"], gridTimes))
Re_gridded_hpd <- t(getHPDMedian(Re_gridded))

times <- mostrecent_decimal - gridTimes*0.1
plotSkyline(times, Re_gridded_hpd, xlab="Date", ylab="Re", type="smooth")   
plotSkyline(times, Re_gridded_hpd, xlab="Date", ylab="Re", type="lines")  

setwd(input_dir)
write.csv(Re_gridded_hpd, file = paste(file_name,".csv",sep = ""), row.names = TRUE)

