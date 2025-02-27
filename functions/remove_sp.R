remove_sp <- function(com, n_sp = 0){
  com_oc <- decostand(com, method = "pa")
  com <- com[,colSums(na.omit(com_oc)) > n_sp]
  return(com)
}