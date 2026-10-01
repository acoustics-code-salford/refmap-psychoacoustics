

## Assign nearest numbers from list --------
numberSnap <- function(numbersIn, numberList, rangeVal){
  
  numberList<-sort(numberList)          # need them sorted
  rangeVal <- rangeVal*1.000001                             # avoid rounding issues
  nearest <- findInterval(numbersIn, numberList - rangeVal) # index of nearest
  nearest <- c(-Inf, numberList)[nearest + 1]  # value of nearest
  diff <- numbersIn - nearest                           # compute errors
  snap <- diff <= rangeVal                               # only snap near numbers
  
  numbersOut <- numbersIn
  numbersOut[snap] <- nearest[snap]                      # snap values to nearest
  
  return(numbersOut)
  
}