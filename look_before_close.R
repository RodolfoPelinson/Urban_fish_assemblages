



#install.packages("renv")
#renv::init()

#Snapshot?
renv::snapshot()
renv::restore()



#Packages
install.packages("usethis")


install.packages("devtools")
library(devtools)
#install.packages("Rcpp")
#devtools::install_github("JenniNiku/gllvm")
install.packages("vegan")
install.packages("glmmTMB", type = "source")
install.packages("emmeans")
install.packages("car")
install.packages("DHARMa")
install.packages("bbmle")
install.packages("TMB")
install.packages("Matrix")
install.packages("randomForest")
install.packages("mgcv")
install.packages("mvabund")

install.packages("DBI")
install.packages("hms")



#install.packages("spaa")
#install.packages("doParallel")
#install.packages("pbapply")


usethis::use_git()
usethis::use_github()
usethis::use_readme_rmd()
