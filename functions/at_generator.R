at_generator <- function(first, spacing, n){
  axis <- first
  for(i in 2:(n)){
    axis[i] <- first + spacing*(i-1)
  }
  return(axis)
}
