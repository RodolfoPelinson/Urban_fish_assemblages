forward_selection_vegan <- function(Y,X,permutations = 9999){
  
  New_X <- X[is.na(rowSums(X)) == FALSE,]
  New_Y <- Y[is.na(rowSums(X)) == FALSE,]
  
  New_X <- data.frame(New_X)
  colnames(New_X) <- colnames(X)
  
  X_rda <- rda(New_Y,New_X)
  p <- anova.cca(X_rda, permutations = permutations)
  if(na.omit(p$`Pr(>F)`) <= 0.05){
    X.R2 <- RsquareAdj(X_rda)$adj.r.squared
    
    res <- try(X_forw <- forward.sel(New_Y, as.matrix(New_X), adjR2thresh = X.R2, nperm = permutations))
    if(inherits(res, "try-error")){
      message("No variables selected")
      New_X_2 <- data.frame(X)
    }else{
      New_X_2 <- data.frame(X[,match(X_forw$variables, colnames(X))])
      colnames(New_X_2) <- colnames(X)[1:nrow(X_forw)]
    }
    
  }
  else{
    message("Forward selection NOT performed. p > 0.05")
    New_X_2 <- data.frame(X)
    X_forw <- p
  }
  if(is.null(New_X_2)){New_X_2 <- data.frame(X)}
  result <- list(forward_results = X_forw,selected_variables = New_X_2)
  return(result)
}
