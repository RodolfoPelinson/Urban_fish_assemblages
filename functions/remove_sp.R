remove_sp <- function(com, n_sp = 0, less_equal = FALSE){
  com_oc <- decostand(com, method = "pa")
  com_more <- com[,colSums(na.omit(com_oc)) > n_sp]
  com_less <- com[,colSums(na.omit(com_oc)) <= n_sp]
  
  if(isTRUE(less_equal)){
    return(com_less)
  }else{
    return(com_more)
  }
}
