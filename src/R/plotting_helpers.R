# plotting_helpers.R
require(ggplot2)

# Rescale ggplot tick marks for variables --------------------------------------

rescaleTicks <- function(dataVar, dataVarScl, sep=3, ax='x'){
  
  if (ax == 'x'){
    outAx <- scale_x_continuous(breaks=seq(floor(min(dataVar)),
                                           ceiling(max(dataVar)),
                                           by=sep) -
                                  (max(dataVar) - max(dataVarScl)),
                                labels=seq(floor(min(dataVar)),
                                           ceiling(max(dataVar)),
                                           by=sep))
    
  } else if (ax == 'y'){
    outAx <- scale_y_continuous(breaks=seq(floor(min(dataVar)),
                                           ceiling(max(dataVar)),
                                           by=sep) -
                                  (max(dataVar) - max(dataVarScl)),
                                labels=seq(floor(min(dataVar)),
                                           ceiling(max(dataVar)),
                                           by=sep))
  } else {
    stop("Invalid axis argument")
  }
  
  return(outAx) 
}
