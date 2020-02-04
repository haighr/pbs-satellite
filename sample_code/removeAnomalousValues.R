removeAnomalousValues <- function(ncdfData, zlim) {
  if (length(zlim) != 2) {
      ## (still) not happy with this warning message
      stop("zlim must be a range of two values; use NA to omit a minimum/maximum")
  }

  ncdfData <- lapply(ncdfData, function(slice, zlim) {
      if (zlim[1] != NA){
          slice[slice < zlim[1]] <- NA
      }
      if (zlim[2] != NA){
          slice[slice > zlim[2]] <- NA
      }

      return(slice)
  }, zlim)

  return(ncdfData)
}
