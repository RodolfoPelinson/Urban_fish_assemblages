text_repositioning <- function(loadings, amount = 0.1) {
 
  x <- loadings[, 1]
  y <- loadings[, 2]
  
  comprimentos <- sqrt(x^2 + y^2)

  
  novo_x <- (((comprimentos+amount)/comprimentos) * x)
  novo_y <- (((comprimentos+amount)/comprimentos) * y)
  
  return(data.frame(x = novo_x, y = novo_y))
}
