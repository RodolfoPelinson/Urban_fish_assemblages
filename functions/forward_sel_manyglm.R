forward_sel_manyglm <- function(y, x, R2more = 0.001, alpha = 0.05, nBoot=999){
  
  R2 <- rep(NA, ncol(x))
  #p <- rep(NA, ncol(x))
  
  y <- mvabund(y)
  
  for(i in 1:ncol(x)){
    
    right_fm <- paste(colnames(x)[i], collapse = " + ")
    formula <- formula(paste("y ~", right_fm))
    
    model <- manyglm(formula = formula, data = x, family = "negative.binomial")
    model_null <- manyglm(formula = y ~ 1, data = x, family = "negative.binomial")
      
    R2[i] <- R2_manyglm(model, model_null)$com_R2$R2
    #anova <- anova.manyglm(model_null, model, nBoot = nBoot)
    #p[i] <- anova$table$`Pr(>Dev)`[2]
  }
  
  names(R2) <- colnames(x)
  #names(p) <- colnames(x)
  
  reorganized_x <- data.frame(x[,order(R2, decreasing = TRUE)])
  colnames(reorganized_x) <- colnames(x)[order(R2, decreasing = TRUE)]
  
  formula_list <- list()
  formula_list[[1]] <- formula("y ~ 1")
  
  right_fm <- paste(colnames(reorganized_x)[1])
  formula <- formula(paste("y ~", right_fm))
  
  formula_list[[2]] <- formula
  
  if(ncol(x) > 1){
    for( i in 2:(ncol(reorganized_x))){
      right_fm <- paste(right_fm, colnames(reorganized_x)[i], sep = " + ")
      formula_list[[i+1]] <- formula(paste("y ~", right_fm))
    }
  }
  
  

  
  result_list <- list()
  
  for( i in 1:ncol(reorganized_x)){
    model <- manyglm(formula = formula_list[[i+1]], data = reorganized_x, family = "negative.binomial")
    model_null <- manyglm(formula = formula_list[[i]], data = reorganized_x, family = "negative.binomial")
    
    R2 <- R2_manyglm(model, model_null)$com_R2$R2
    anova <- anova.manyglm(model_null, model, nBoot = nBoot)
    p <- anova$table$`Pr(>Dev)`[2]
    df.diff <- anova$table$Df.diff[2]
    Dev <- anova$table$Dev[2]
    
    result_list[[i]]<-c(df.diff = df.diff, Dev = Dev, R2 = R2, p = p)
    names(result_list)[i] <- colnames(reorganized_x)[i]
    
    if(R2 < R2more | p > alpha){
      break
    }
    
  }

result <- t(as.data.frame(result_list))

new_x_id <- match(rownames(result), colnames(x))

if(ncol(x)>1){
  new_x <- data.frame(x[,new_x_id[-length(new_x_id)]])
  colnames(new_x) <- colnames(x)[new_x_id[-length(new_x_id)]]
}else{
  new_x <- data.frame(x[,new_x_id])
  colnames(new_x) <- colnames(x)[new_x_id]
}


result_list <- list(result = result, new_x = new_x)

print(result)  

return(result_list)

}


