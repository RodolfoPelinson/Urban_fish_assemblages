#This function compares models with and without quadratic terms when quad is true


forward_sel_manyglm2 <- function(y, x, R2more = 0.001, alpha = 0.05, nBoot=999, quad = FALSE){
  
  R2 <- rep(NA, ncol(x))
  #p <- rep(NA, ncol(x))
  
  y <- mvabund(y)
  
  new_x <- x
  new_result <- NULL
  
  if(isTRUE(quad)){
    
    message("testing for quadratic effects...")
    
    
    x_quad <- x^2
    colnames(x_quad) <- paste(colnames(x_quad), "squared", sep = "_")
    
    for(i in 1:ncol(x)){
      
      right_fm <- paste(c(colnames(x)[i], colnames(x_quad)[i]), collapse = " + ")
      right_fm_null <- paste(c(colnames(x)[i]), collapse = " + ")
      formula <- formula(paste("y ~", right_fm))
      formula_null <- formula(paste("y ~", right_fm_null))
      
      model <- manyglm(formula = formula, data = data.frame(x, x_quad), family = "negative.binomial")
      model_null <- manyglm(formula = formula_null, data = x, family = "negative.binomial")
      
      R2[i] <- R2_manyglm(model, model_null)$com_R2$R2
      #anova <- anova.manyglm(model_null, model, nBoot = nBoot)
      #p[i] <- anova$table$`Pr(>Dev)`[2]
      
      
    }
    names(R2) <- colnames(x)
    
    
    reorganized_x <- data.frame(x[,order(R2, decreasing = TRUE)])
    colnames(reorganized_x) <- colnames(x)[order(R2, decreasing = TRUE)]
    
    reorganized_x_quad <- data.frame(x_quad[,order(R2, decreasing = TRUE)])
    colnames(reorganized_x_quad) <- colnames(x_quad)[order(R2, decreasing = TRUE)]
    
    formula_list <- list()
    formula_list_null <- list()
    
    right_fm <- paste(c(colnames(reorganized_x[1]), colnames(reorganized_x_quad[1])), collapse = " + ")
    formula <- formula(paste("y ~", right_fm))
    
    right_fm_null <- paste(c(colnames(reorganized_x[1])), collapse = " + ")
    formula_null <- formula(paste("y ~", right_fm_null))
    
    formula_list[[1]] <- formula
    formula_list_null[[1]] <- formula_null
    
    
    if(ncol(x) > 1){
      for( i in 2:(ncol(reorganized_x))){
        new_right_fm <- paste(c(colnames(reorganized_x[i]), colnames(reorganized_x_quad[i])), collapse = " + ")
        right_fm <- paste(right_fm, new_right_fm, sep = " + ")
        formula_list[[i]] <- formula(paste("y ~", right_fm))
        
        new_right_fm_null <- paste(c(colnames(reorganized_x[i])), collapse = " + ")
        right_fm_null <- paste(c(colnames(reorganized_x[i-1]), colnames(reorganized_x_quad[i-1]), new_right_fm_null), collapse = " + ")
        formula_list_null[[i]] <- formula(paste("y ~", right_fm_null))
        
      }
    }
    
    
    
    data_x <- data.frame(reorganized_x, reorganized_x_quad)
    
    
    
    result_list <- list()
    
    for( i in 1:ncol(reorganized_x)){
      model <- manyglm(formula = formula_list[[i]], data = data_x, family = "negative.binomial")
      model_null <- manyglm(formula = formula_list_null[[i]], data = data_x, family = "negative.binomial")
      
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
    
    print(result)
    
    new_x_id <- match(rownames(result), colnames(x))
    
    if(ncol(x)>1){
      new_x <- data.frame(x[,new_x_id[-length(new_x_id)]])
      colnames(new_x) <- colnames(x)[new_x_id[-length(new_x_id)]]
    }else{
      new_x <- data.frame(x[,new_x_id])
      colnames(new_x) <- colnames(x)[new_x_id]
    }
    
    
    
    new_x_id_quad <- match(rownames(result), colnames(x))
    
    if(ncol(x_quad)>1){
      new_x_quad <- data.frame(x_quad[,new_x_id_quad[-length(new_x_id_quad)]])
      colnames(new_x_quad) <- colnames(x_quad)[new_x_id_quad[-length(new_x_id_quad)]]
    }else{
      new_x_quad <- data.frame(x_quad[,new_x_id_quad])
      colnames(new_x_quad) <- colnames(x_quad)[new_x_id_quad]
    }
    
    new_x <- data.frame(new_x, new_x_quad)
    new_x <- new_x[,order(colnames(new_x))]
    
    if(ncol(new_x) == 0){
      message("No quadratic effects were found, testing for linear ones...")
    }else{
      message("testing for linear effects...")
    }
    
    
    
    
    if(ncol(new_x) > 0){
      
      new_reorganized_x <- reorganized_x[(nrow(result)):ncol(reorganized_x)]
      
      #Reorganizar as variaveis que sobraram #################3
      base_formula <- formula_list[nrow(result)-1]
      
      for(i in 1:ncol(new_reorganized_x)){
        
        right_fm <- paste(colnames(new_reorganized_x)[i], collapse = " + ")
        formula <- formula(paste(paste(base_formula), right_fm, sep =" + "))
        
        model <- manyglm(formula = formula, data = data.frame(reorganized_x, reorganized_x_quad), family = "negative.binomial")
        model_null <- manyglm(formula = base_formula[[1]], data = data.frame(reorganized_x, reorganized_x_quad), family = "negative.binomial")
        
        R2[i] <- R2_manyglm(model, model_null)$com_R2$R2
        #anova <- anova.manyglm(model_null, model, nBoot = nBoot)
        #p[i] <- anova$table$`Pr(>Dev)`[2]
      }
      
      names(R2) <- colnames(new_reorganized_x)
      #names(p) <- colnames(x)
      
      new_new_reorganized_x <- data.frame(new_reorganized_x[,order(R2, decreasing = TRUE)])
      colnames(new_new_reorganized_x) <- colnames(new_reorganized_x)[order(R2, decreasing = TRUE)]
      
      new_reorganized_x <- new_new_reorganized_x
      
      #########################################################
      
      
      
      # base_formula <- formula_list[nrow(result)]
      
      
      formula_list <- list()
      formula_list[[1]] <- base_formula[[1]]
      
      right_fm <- paste(colnames(new_reorganized_x)[1])
      formula <- formula(paste(paste(base_formula), right_fm, sep =" + "))
      
      
      
      
      
      
      formula_list[[2]] <- formula
      
      if(ncol(new_reorganized_x) > 1){
        for( i in 2:(ncol(new_reorganized_x))){
          right_fm <- paste(right_fm, colnames(new_reorganized_x)[i], sep = " + ")
          formula_list[[i+1]] <- formula(paste(base_formula, right_fm, sep = " + "))
        }
      }
      
      
      
      data_x <- data.frame(new_reorganized_x, new_x)
      
      
      result_list <- list()
      
      for( i in 1:ncol(new_reorganized_x)){
        model <- manyglm(formula = formula_list[[i+1]], data = data_x, family = "negative.binomial")
        model_null <- manyglm(formula = formula_list[[i]], data = data_x, family = "negative.binomial")
        
        R2 <- R2_manyglm(model, model_null)$com_R2$R2
        anova <- anova.manyglm(model_null, model, nBoot = nBoot)
        p <- anova$table$`Pr(>Dev)`[2]
        df.diff <- anova$table$Df.diff[2]
        Dev <- anova$table$Dev[2]
        
        result_list[[i]]<-c(df.diff = df.diff, Dev = Dev, R2 = R2, p = p)
        names(result_list)[i] <- colnames(new_reorganized_x)[i]
        
        if(R2 < R2more | p > alpha){
          break
        }
        
      }
      
      new_result <- t(as.data.frame(result_list))
      
      print(new_result)
      
      new_new_x_id <- match(rownames(new_result), colnames(x))
      
      if(ncol(x)>1){
        new_new_x <- data.frame(x[,new_new_x_id[-length(new_new_x_id)]])
        colnames(new_new_x) <- colnames(x)[new_new_x_id[-length(new_new_x_id)]]
      }else{
        new_new_x <- data.frame(x[,new_new_x_id])
        colnames(new_new_x) <- colnames(x)[new_new_x_id]
      }
      
      if(ncol(new_new_x) == 0){
        message("No linear effects were found")
      }
      
      new_x<-data.frame(new_x, new_new_x)
      
    }
    
    
    
    
    
  }
  
  
  
  ########################################################################
  
  
  
  
  
  if(quad == FALSE | ncol(new_x) == 0){
    
    
    
    
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
    
    
    
    data_x <- data.frame(reorganized_x)
    
    
    result_list <- list()
    
    for( i in 1:ncol(reorganized_x)){
      model <- manyglm(formula = formula_list[[i+1]], data = data_x, family = "negative.binomial")
      model_null <- manyglm(formula = formula_list[[i]], data = data_x, family = "negative.binomial")
      
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
    
    print(result)
    
    new_x_id <- match(rownames(result), colnames(x))
    
    if(ncol(x)>1){
      new_x <- data.frame(x[,new_x_id[-length(new_x_id)]])
      colnames(new_x) <- colnames(x)[new_x_id[-length(new_x_id)]]
    }else{
      new_x <- data.frame(x[,new_x_id])
      colnames(new_x) <- colnames(x)[new_x_id]
    }
    
    if(ncol(new_x) == 0){
      message("No linear effects were found, adding the first predictor to new_x")
      new_x<-reorganized_x[1]
      
    }
    
    #  new_x<-reorganized_x[1]
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  if(is.null(new_result)){
    result_list <- list(result_linear = result, new_x = new_x)
  }else{
    result_list <- list(result_quad = result, result_linear = new_result, new_x = new_x)
    
  }
  
  return(result_list)
  
}


