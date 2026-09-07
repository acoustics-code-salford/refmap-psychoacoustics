# plotting_helpers.R
require(ggplot2)

# Rescale ggplot tick marks for variables --------------------------------------

rescaleTicks <- function(dataVar, dataVarScl, sep=3, ax='x', sd.scale=FALSE){

  data_mean <- mean(dataVar, na.rm=TRUE)
  scale_factor <- if (sd.scale) sd(dataVar, na.rm=TRUE) else 1
  
  dataVar_min <- min(dataVar, na.rm=TRUE)
  dataVar_max <- max(dataVar, na.rm=TRUE)
  
  ticks <- seq(floor(dataVar_min), ceiling(dataVar_max), by=sep)
  
  breaks <- (ticks - data_mean) / scale_factor
  
  if (ax == 'x'){
    outAx <- scale_x_continuous(breaks=breaks, labels=ticks)
  } else if (ax == 'y'){
    outAx <- scale_y_continuous(breaks=breaks, labels=ticks)
  } else {
    stop("Invalid axis argument")
  }
  
  return(outAx) 
}
