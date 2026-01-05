



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

install.packages("adespatial")
install.packages("ade4")

#install.packages("spaa")
#install.packages("doParallel")
#install.packages("pbapply")


install.packages("VennDiagram")

VennDiagram

install.packages("usethis")

#Para configurar git pela primeira vez no PC
#usethis::use_git_config(
#  user.name = "Rodolfo Pelinson",
#  user.email = "rodolfopelinson@gmail.com"
#)

usethis::git_sitrep()


usethis::use_git()

#Para conectar o repositorio a um já existente no github
usethis::git_remotes_set(
  name = "origin",
  url  = "https://github.com/RodolfoPelinson/Urban_fish_assemblages.git"
)

git remote add origin https://github.com/RodolfoPelinson/Urban_fish_assemblages.git

usethis::use_github()
usethis::use_readme_rmd()
