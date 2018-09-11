# 2_size_spectra

library(plyr)
# function to bin and center data
bin_and_center <- function(data, var, breaks, ...){
  # data is a data frame
  # var is a string, and is the name of a column in data which you want to bin
  # breaks controls the number of bins as defined in hist() 
    # See ?hist for details
  
  # bin values using hist()
  binned_hist = hist(data[[var]], 
                     breaks = breaks, # need to predefine breaks
                     # e.g. Log2 breaks = 2^seq(min, max) 
                     # Log10 breaks = 10^seq(min, max)
                     include.lowest = TRUE, plot = FALSE)
  # calculate "left" and "right" edge of bins
  breaks_orig = binned_hist$breaks[1:(length(breaks)-1)]
  breaks_offset = binned_hist$breaks[2:length(breaks)]
  # total bin width = right edge - left edge
  break_width = breaks_offset - breaks_orig
  count = binned_hist$counts 
  dataout = data.frame(
    # normalize counts =count / width (White et al 1997)
    log_count_corrected = log10(count / break_width),
    # original midpoint of bin log10 transformed
    log_mids = log10(binned_hist$mids),
    log_mids_center = NA)
  # remove bins with 0 counts
  # -Inf comes from log10(count / break_width) above
  dataout = dataout[dataout$log_count_corrected !=-Inf,]
  # recenter data at x=0
    mid_row = ceiling(nrow(dataout)/2)
    # subtract value of mid row from all mids
    dataout$log_mids_center = 
      dataout[,"log_mids"] - dataout[mid_row,"log_mids"]
  dataout
}

# data in ####
# read in estimated dry weights
# dw is in units of grams
dataset <- read.csv("data/estimated_dw.csv",
                    stringsAsFactors = FALSE)
# read in pca.axis
pca.axis <- read.csv("results/pca1_axis.csv",
                     stringsAsFactors = FALSE)

# bin dataset ####
# predefine log2 width break points 
breaks <- 2^seq(-22,-2)

# make sure that range of breaks is
# greater than range of data
test1 <- min(breaks) < min(dataset$dw) & max(breaks) > max(dataset$dw)

if(test1 == TRUE){
  # applies bin_and_center() to the dataset, after subsetting by "site" variable
    binned <- ddply(dataset, .(site),
                bin_and_center,
                var = "dw",
                breaks = breaks)
}else{
  "Breaks does not cover range of data"
}

# adding pca axes
binned <- merge(pca.axis, binned, by = "site")

# save results
write.csv(binned, "results/binned_size_spectra.csv",
          row.names = FALSE)
