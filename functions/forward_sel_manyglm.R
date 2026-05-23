#This function compares models with quadratic terms to models with no terms

forward_sel_manyglm <- function(y, x, alpha = 0.05, nBoot=999, quad = FALSE, adj_R2 = TRUE){
  
  R2 <- rep(NA, ncol(x))
  p <- rep(NA, ncol(x))
  
  y <- mvabund(y)
  
  new_x <- x
  new_result <- NULL
  
  if(isTRUE(quad)){
    
    message("testing for quadratic effects...")
    message("testing for Global Model...")
    x_quad <- x^2
    colnames(x_quad) <- paste(colnames(x_quad), "squared", sep = "_")
    
    right_fm <- paste(c(colnames(x), colnames(x_quad)), collapse = " + ")
    formula <- formula(paste("y ~", right_fm))
    global_model <- manyglm(formula = formula, data = data.frame(x, x_quad), family = "negative.binomial")
    model_null <- manyglm(formula = y ~ 1, data = x, family = "negative.binomial")
    anova_global <- anova.manyglm(model_null, global_model, nBoot = nBoot)
    p_global <- anova_global$table$`Pr(>Dev)`[2]
    
    if(isTRUE(adj_R2)){
      R2_global <- R2_manyglm(global_model, model_null)$com_R2$adj_R2
    }else{
      R2_global <- R2_manyglm(global_model, model_null)$com_R2$R2
    }
    
    
    if(p_global > alpha){
      message(paste("Global quadratic model is not significant with p value of "), p_global, paste(" and R2 of ", R2_global))
      
      message("Adding all predictor to new_x")
      
      if(isTRUE(quad)){
        new_x<-data.frame(x, x_quad)
        
      }else{
        new_x<-data.frame(x)
      }
      
      result <- list(anova_global = anova_global,
                     R2_global = R2_global)
      
      
      
    }else{
      
      message(paste("Global quadratic model is significant with p value of "), p_global, paste(" and R2 of ", R2_global))
      
      message("Executing forward selection...")
      
      for(i in 1:ncol(x)){
        
        right_fm <- paste(c(colnames(x)[i], colnames(x_quad)[i]), collapse = " + ")
        formula <- formula(paste("y ~", right_fm))
        
        model <- manyglm(formula = formula, data = data.frame(x, x_quad), family = "negative.binomial")
        model_null <- manyglm(formula = y ~ 1, data = x, family = "negative.binomial")
        
        
        if(isTRUE(adj_R2)){
          R2[i] <- R2_manyglm(model, model_null)$com_R2$adj_R2
        }else{
          R2[i] <- R2_manyglm(model, model_null)$com_R2$R2
        }
        #anova <- anova.manyglm(model_null, model, nBoot = nBoot)
        #p[i] <- anova$table$`Pr(>Dev)`[2]
        
        
      }
      names(R2) <- colnames(x)
      #names(p) <- colnames(x)
      
      #significant_x <- data.frame(x[,p <= alpha])
      #significant_x_quad <- data.frame(x_quad[,p <= alpha])
      #significant_R2 <- R2[p <= alpha]
      
      #reorganized_x <- data.frame(significant_x[,order(significant_R2, decreasing = TRUE)])
      #colnames(reorganized_x) <- colnames(x)[p <= alpha][order(significant_R2, decreasing = TRUE)]
      
      #reorganized_x_quad <- data.frame(significant_x_quad[,order(significant_R2, decreasing = TRUE)])
      #colnames(reorganized_x_quad) <- colnames(x_quad)[p <= alpha][order(significant_R2, decreasing = TRUE)]

      reorganized_x <- data.frame(x[,order(R2, decreasing = TRUE)])
      colnames(reorganized_x) <- colnames(x)[order(R2, decreasing = TRUE)]
      
      reorganized_x_quad <- data.frame(x_quad[,order(R2, decreasing = TRUE)])
      colnames(reorganized_x_quad) <- colnames(x_quad)[order(R2, decreasing = TRUE)]
      
      formula_list <- list()
      formula_list[[1]] <- formula("y ~ 1")
      
      right_fm <- paste(c(colnames(reorganized_x[1]), colnames(reorganized_x_quad[1])), collapse = " + ")
      formula <- formula(paste("y ~", right_fm))
      
      formula_list[[2]] <- formula
      
      if(ncol(reorganized_x) > 1){
        for( i in 2:(ncol(reorganized_x))){
          new_right_fm <- paste(c(colnames(reorganized_x[i]), colnames(reorganized_x_quad[i])), collapse = " + ")
          right_fm <- paste(right_fm, new_right_fm, sep = " + ")
          formula_list[[i+1]] <- formula(paste("y ~", right_fm))
        }
      }
      
      
      
      data_x <- data.frame(reorganized_x, reorganized_x_quad)
      
      
      
      
      result_list <- list()
      
      for( i in 1:ncol(reorganized_x)){
        model <- manyglm(formula = formula_list[[i+1]], data = data_x, family = "negative.binomial")
        model_base <- manyglm(formula = formula_list[[i]], data = data_x, family = "negative.binomial")
        
        if(isTRUE(adj_R2)){
          R2 <- R2_manyglm(model, model_null)$com_R2$adj_R2
        }else{
          R2 <- R2_manyglm(model, model_null)$com_R2$R2
        }
        
        anova <- anova.manyglm(model_base, model, nBoot = nBoot)
        p <- anova$table$`Pr(>Dev)`[2]
        df.diff <- anova$table$Df.diff[2]
        Dev <- anova$table$Dev[2]
        
        result_list[[i]]<-c(df.diff = df.diff, Dev = Dev, R2 = R2, p = p)
        names(result_list)[i] <- colnames(reorganized_x)[i]
        
        if(R2 > R2_global | p > alpha){
          break
        }
        
      }
      
      result <- t(as.data.frame(result_list))
      
      print(result)
      
      new_x_id <- match(rownames(result), colnames(x))
      
      
      if(ncol(x)>1 & result[nrow(result),4] < alpha){
        new_x <- data.frame(x[,new_x_id])
        colnames(new_x) <- colnames(x)[new_x_id]
      }
      
      if(ncol(x)>1 & result[nrow(result),4] > alpha){
        new_x <- data.frame(x[,new_x_id[-length(new_x_id)]])
        colnames(new_x) <- colnames(x)[new_x_id[-length(new_x_id)]]
      }
      
      if(ncol(x)<=1){
        new_x <- data.frame(x[,new_x_id])
        colnames(new_x) <- colnames(x)[new_x_id]
      }
      
      
      
      new_x_id_quad <- match(rownames(result), colnames(x))
      
      if(ncol(x_quad)>1 & result[nrow(result),4] < alpha){
        new_x_quad <- data.frame(x_quad[,new_x_id_quad])
        colnames(new_x_quad) <- colnames(x_quad)[new_x_id_quad]
      }
      
      if(ncol(x_quad)>1 & result[nrow(result),4] > alpha){
        new_x_quad <- data.frame(x_quad[,new_x_id_quad[-length(new_x_id_quad)]])
        colnames(new_x_quad) <- colnames(x_quad)[new_x_id_quad[-length(new_x_id_quad)]]
      }
      
      if(ncol(x_quad) <= 1){
        new_x_quad <- data.frame(x_quad[,new_x_id_quad])
        colnames(new_x_quad) <- colnames(x_quad)[new_x_id_quad]
      }
      
      new_x_lin <- new_x
      new_x <- data.frame(new_x, new_x_quad)
      new_x <- new_x[,order(colnames(new_x))]
      
      if(ncol(new_x) == 0){
        
      }else{
        message("testing for linear effects...")
      }
      
      
      if((ncol(new_x)/2) == ncol(x)){
        message("No remaining linear effects to be tested.")
      }else{
        if(ncol(new_x) > 0){
          
          new_reorganized_x <- data.frame(x[,colnames(x) %in% colnames(new_x_lin) == FALSE])
          colnames(new_reorganized_x) <- colnames(x)[colnames(x) %in% colnames(new_x_lin) == FALSE]
          
          #Reorganizar as variaveis que sobraram #################3
          base_formula <- formula_list[nrow(result)]
          
          R2 <- rep(NA, ncol(new_reorganized_x))
          p <- rep(NA, ncol(new_reorganized_x))
          
          for(i in 1:ncol(new_reorganized_x)){
            
            right_fm <- paste(colnames(new_reorganized_x)[i], collapse = " + ")
            formula <- formula(paste(paste(base_formula), right_fm, sep =" + "))
            
            model <- manyglm(formula = formula, data = data.frame(x,x_quad), family = "negative.binomial")
            model_base <- manyglm(formula = base_formula[[1]], data = data.frame(x,x_quad), family = "negative.binomial")
            
            if(isTRUE(adj_R2)){
              R2[i] <- R2_manyglm(model, model_null)$com_R2$adj_R2
            }else{
              R2[i] <- R2_manyglm(model, model_null)$com_R2$R2
            }
            
            #anova <- anova.manyglm(model_base, model, nBoot = nBoot)
            #p[i] <- anova$table$`Pr(>Dev)`[2]
          }
          
          names(R2) <- colnames(new_reorganized_x)
          #names(p) <- colnames(new_reorganized_x)
          
          #significant_new_reorganized_x <- data.frame(new_reorganized_x[,p <= alpha])
          #significant_R2 <- R2[p <= alpha]
          if(ncol(new_reorganized_x)>0){
            
          #if(ncol(significant_new_reorganized_x)>0){
            #significant_new_reorganized_x <- data.frame(significant_new_reorganized_x[,order(significant_R2, decreasing = TRUE)])
            #colnames(significant_new_reorganized_x) <- colnames(new_reorganized_x)[p <= alpha][order(significant_R2, decreasing = TRUE)]
            
            #new_new_reorganized_x <- data.frame(significant_new_reorganized_x[,order(significant_R2, decreasing = TRUE)])
            #colnames(new_new_reorganized_x) <- colnames(significant_new_reorganized_x)[order(significant_R2, decreasing = TRUE)]
            
            new_new_reorganized_x <- data.frame(new_reorganized_x[,order(R2, decreasing = TRUE)])
            colnames(new_new_reorganized_x) <- colnames(new_reorganized_x)[order(R2, decreasing = TRUE)]
            
            
            new_reorganized_x <- new_new_reorganized_x
            
            
            
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
              
              if(isTRUE(adj_R2)){
                R2 <- R2_manyglm(model, model_null)$com_R2$adj_R2
              }else{
                R2 <- R2_manyglm(model, model_null)$com_R2$R2
              }
              anova <- anova.manyglm(model_null, model, nBoot = nBoot)
              p <- anova$table$`Pr(>Dev)`[2]
              df.diff <- anova$table$Df.diff[2]
              Dev <- anova$table$Dev[2]
              
              result_list[[i]]<-c(df.diff = df.diff, Dev = Dev, R2 = R2, p = p)
              names(result_list)[i] <- colnames(new_reorganized_x)[i]
              
              if(R2 < R2_global | p > alpha){
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
            
            
    
          }else{
            message("No linear effects were found")
          }
          

          
          #########################################################
          
          
          

          
        }
      }
    }

    

    

    
    
    
    
    
  }
  
  
  
  ########################################################################
  
  
  #if(quad == FALSE | ncol(new_x) == 0){
  #if(quad == TRUE & ncol(new_x) == 0){
    
  
  
  if(quad == FALSE){
    
      message("testing for linear effects...")
      message("testing for Global Model...")
      
      right_fm <- paste(colnames(x), collapse = " + ")
      formula <- formula(paste("y ~", right_fm))
      global_model <- manyglm(formula = formula, data = data.frame(x), family = "negative.binomial")
      model_null <- manyglm(formula = y ~ 1, data = x, family = "negative.binomial")
      anova_global <- anova.manyglm(model_null, global_model, nBoot = nBoot)
      p_global <- anova_global$table$`Pr(>Dev)`[2]
      
      if(isTRUE(adj_R2)){
        R2_global <- R2_manyglm(global_model, model_null)$com_R2$adj_R2
      }else{
        R2_global <- R2_manyglm(global_model, model_null)$com_R2$R2
      }
      
      
    if(p_global > alpha){
      message(paste("Global linear model is not significant with p value of "), p_global, paste(" and R2 of ", R2_global))
      
      message("Adding all predictor to new_x")
      
      new_x<-data.frame(x)
      
      
      result <- list(anova_global = anova_global,
                     R2_global = R2_global)
    
      
    }else{
      
      message(paste("Global linear model is significant with p value of "), p_global, paste(" and R2 of ", R2_global))
      
      message("Executing forward selection...")
      
      for(i in 1:ncol(x)){
        
        right_fm <- paste(colnames(x)[i], collapse = " + ")
        formula <- formula(paste("y ~", right_fm))
        
        #right_fm_global <- paste(colnames(x), collapse = " + ")
        #formula_global <- formula(paste("y ~", right_fm_global))
        
        #model_global <- manyglm(formula = formula_global, data = data.frame(x), family = "negative.binomial")
        
        
        model <- manyglm(formula = formula, data = data.frame(x), family = "negative.binomial")
        model_null <- manyglm(formula = y ~ 1, data = x, family = "negative.binomial")
        
        
        
        if(isTRUE(adj_R2)){
          R2[i] <- R2_manyglm(model, model_null)$com_R2$adj_R2
        }else{
          R2[i] <- R2_manyglm(model, model_null)$com_R2$R2
        }
        #anova <- anova.manyglm(model_conditional, model_global, nBoot = nBoot)
        #p[i] <- anova$table$`Pr(>Dev)`[2]
      }
      
      names(R2) <- colnames(x)
      #names(p) <- colnames(x)
      
      
      #significant_x <- data.frame(x[,p <= alpha])
      #significant_R2 <- R2[p <= alpha]
      
      #reorganized_x <- data.frame(significant_x[, order(significant_R2, decreasing = TRUE)])
      #colnames(reorganized_x) <- colnames(x)[p <= alpha][order(significant_R2, decreasing = TRUE)]
      
      reorganized_x <- data.frame(x[, order(R2, decreasing = TRUE)])
      colnames(reorganized_x) <- colnames(x)[order(R2, decreasing = TRUE)]
      
      formula_list <- list()
      formula_list[[1]] <- formula("y ~ 1")
      
      right_fm <- paste(colnames(reorganized_x)[1])
      formula <- formula(paste("y ~", right_fm))
      
      formula_list[[2]] <- formula
      
      if(ncol(reorganized_x) > 1){
        for( i in 2:(ncol(reorganized_x))){
          right_fm <- paste(right_fm, colnames(reorganized_x)[i], sep = " + ")
          formula_list[[i+1]] <- formula(paste("y ~", right_fm))
        }
      }
      
      
      
      data_x <- data.frame(reorganized_x)
      
      
      result_list <- list()
      
      for( i in 1:ncol(reorganized_x)){
        model <- manyglm(formula = formula_list[[i+1]], data = data_x, family = "negative.binomial")
        model_base <- manyglm(formula = formula_list[[i]], data = data_x, family = "negative.binomial")
        
        if(isTRUE(adj_R2)){
          R2 <- R2_manyglm(model, model_null)$com_R2$adj_R2
        }else{
          R2 <- R2_manyglm(model, model_null)$com_R2$R2
        }
        anova <- anova.manyglm(model_base, model, nBoot = nBoot)
        p <- anova$table$`Pr(>Dev)`[2]
        df.diff <- anova$table$Df.diff[2]
        Dev <- anova$table$Dev[2]
        
        result_list[[i]]<-c(df.diff = df.diff, Dev = Dev, R2 = R2, p = p)
        names(result_list)[i] <- colnames(reorganized_x)[i]
        
        if(R2 > R2_global | p > alpha){
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
        message("Only the global model is significant, adding all predictor to new_x")
        
        if(isTRUE(quad)){
          new_x<-data.frame(x, x_quad)
          
        }else{
          new_x<-data.frame(x)
        }
      }
    }
    
    
    

    

    
    #  new_x<-reorganized_x[1]
    
  }
  
  
  
  if(quad == TRUE & ncol(new_x) == 0){
  
    message("testing for linear effects...")

    
      for(i in 1:ncol(x)){
        
        right_fm <- paste(colnames(x)[i], collapse = " + ")
        formula <- formula(paste("y ~", right_fm))
        
        #right_fm_global <- paste(colnames(x), collapse = " + ")
        #formula_global <- formula(paste("y ~", right_fm_global))
        
        #model_global <- manyglm(formula = formula_global, data = data.frame(x), family = "negative.binomial")
        
        
        model <- manyglm(formula = formula, data = data.frame(x), family = "negative.binomial")
        model_null <- manyglm(formula = y ~ 1, data = x, family = "negative.binomial")
        
        
        
        if(isTRUE(adj_R2)){
          R2[i] <- R2_manyglm(model, model_null)$com_R2$adj_R2
        }else{
          R2[i] <- R2_manyglm(model, model_null)$com_R2$R2
        }
        #anova <- anova.manyglm(model_conditional, model_global, nBoot = nBoot)
        #p[i] <- anova$table$`Pr(>Dev)`[2]
      }
      
      names(R2) <- colnames(x)
      #names(p) <- colnames(x)
      
      
      #significant_x <- data.frame(x[,p <= alpha])
      #significant_R2 <- R2[p <= alpha]
      
      #reorganized_x <- data.frame(significant_x[, order(significant_R2, decreasing = TRUE)])
      #colnames(reorganized_x) <- colnames(x)[p <= alpha][order(significant_R2, decreasing = TRUE)]
      
      reorganized_x <- data.frame(x[, order(R2, decreasing = TRUE)])
      colnames(reorganized_x) <- colnames(x)[order(R2, decreasing = TRUE)]
      
      formula_list <- list()
      formula_list[[1]] <- formula("y ~ 1")
      
      right_fm <- paste(colnames(reorganized_x)[1])
      formula <- formula(paste("y ~", right_fm))
      
      formula_list[[2]] <- formula
      
      if(ncol(reorganized_x) > 1){
        for( i in 2:(ncol(reorganized_x))){
          right_fm <- paste(right_fm, colnames(reorganized_x)[i], sep = " + ")
          formula_list[[i+1]] <- formula(paste("y ~", right_fm))
        }
      }
      
      
      
      data_x <- data.frame(reorganized_x)
      
      
      result_list <- list()
      
      for( i in 1:ncol(reorganized_x)){
        model <- manyglm(formula = formula_list[[i+1]], data = data_x, family = "negative.binomial")
        model_base <- manyglm(formula = formula_list[[i]], data = data_x, family = "negative.binomial")
        
        if(isTRUE(adj_R2)){
          R2 <- R2_manyglm(model, model_null)$com_R2$adj_R2
        }else{
          R2 <- R2_manyglm(model, model_null)$com_R2$R2
        }
        anova <- anova.manyglm(model_base, model, nBoot = nBoot)
        p <- anova$table$`Pr(>Dev)`[2]
        df.diff <- anova$table$Df.diff[2]
        Dev <- anova$table$Dev[2]
        
        result_list[[i]]<-c(df.diff = df.diff, Dev = Dev, R2 = R2, p = p)
        names(result_list)[i] <- colnames(reorganized_x)[i]
        
        if(R2 > R2_global | p > alpha){
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
        message("Only the global model is significant, adding all predictor to new_x")
        
        if(isTRUE(quad)){
          new_x<-data.frame(x, x_quad)
          
        }else{
          new_x<-data.frame(x)
        }
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