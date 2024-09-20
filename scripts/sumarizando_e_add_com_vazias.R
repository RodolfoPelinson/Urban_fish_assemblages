abundance_orig <- read.csv("data/abundancia_peixes_por_trecho.csv")
delineamento <- read.csv("data/delineamento.csv")


abundancia <- subset(abundance_orig, select = -c(Poecilia.cf.))

#Adicionando riachos sem sem peixe
bacia_sem_peixe <- setdiff(delineamento$bacia_id, abundancia$Riacho)
bacias_sem_peixe <- matrix(0, ncol = ncol(abundancia), nrow = length(bacia_sem_peixe)*5)
colnames(bacias_sem_peixe) <- colnames(abundancia)


Trechos <- c(1:5)
Riachos <- rep(bacia_sem_peixe[1],5)

for(i in 2:length(bacia_sem_peixe)){
  Riacho <- rep(bacia_sem_peixe[i],5)
  Trecho <- c(1:5)
  Riachos <- c(Riachos,Riacho)
  Trechos <- c(Trechos,Trecho)
}

bacias_sem_peixe <- as.data.frame(bacias_sem_peixe)

bacias_sem_peixe$Riacho <- Riachos
bacias_sem_peixe$Trecho <- Trechos

abundancia_trechos <- rbind(abundancia, bacias_sem_peixe)

com <- abundancia_trechos


################################## sumarizando
Bacias <- unique(abundancia_trechos$Riacho)
comunidade_bacia <- abundancia_trechos[abundancia_trechos == Bacias[1],-c(1:2)]
comunidade <- apply(comunidade_bacia, 2, sum)
abundancia_bacia <- c(Riacho = Bacias[1], comunidade)

for(i in 2:length(Bacias)){
  comunidade_bacia <- abundancia_trechos[abundancia_trechos == Bacias[i],-c(1:2)]
  comunidade <- apply(comunidade_bacia, 2, sum)
  abundancia_bacia <- rbind(abundancia_bacia, c(Bacias[i], comunidade))
}

rownames(abundancia_bacia) <- abundancia_bacia[,1]
abundancia_bacia <- as.data.frame(abundancia_bacia[,-1])

abundancia_bacia <- apply(abundancia_bacia, 2, as.numeric); rownames(abundancia_bacia) <- Bacias
abundancia_bacia <- as.data.frame(abundancia_bacia)

com_sum <- abundancia_bacia

########################################




#################################################
############### BIOMASSA ########################
#################################################

biomass_orig <- read.csv("data/biomassa_peixes_por_trecho.csv")
delineamento <- read.csv("data/delineamento.csv")


biomassa <- subset(biomass_orig, select = -c(Poecilia.cf.))

#Adicionando riachos sem sem peixe
bacia_sem_peixe <- setdiff(delineamento$bacia_id, biomassa$Riacho)
bacias_sem_peixe <- matrix(0, ncol = ncol(biomassa), nrow = length(bacia_sem_peixe)*5)
colnames(bacias_sem_peixe) <- colnames(biomassa)


Trechos <- c(1:5)
Riachos <- rep(bacia_sem_peixe[1],5)

for(i in 2:length(bacia_sem_peixe)){
  Riacho <- rep(bacia_sem_peixe[i],5)
  Trecho <- c(1:5)
  Riachos <- c(Riachos,Riacho)
  Trechos <- c(Trechos,Trecho)
}

bacias_sem_peixe <- as.data.frame(bacias_sem_peixe)

bacias_sem_peixe$Riacho <- Riachos
bacias_sem_peixe$Trecho <- Trechos

biomassa_trechos <- rbind(biomassa, bacias_sem_peixe)

com_biomass <- biomassa_trechos


################################## sumarizando
Bacias <- unique(biomassa_trechos$Riacho)
comunidade_bacia <- biomassa_trechos[biomassa_trechos == Bacias[1],-c(1:2)]
comunidade <- apply(comunidade_bacia, 2, sum)
biomassa_bacia <- c(Riacho = Bacias[1], comunidade)

for(i in 2:length(Bacias)){
  comunidade_bacia <- biomassa_trechos[biomassa_trechos == Bacias[i],-c(1:2)]
  comunidade <- apply(comunidade_bacia, 2, sum)
  biomassa_bacia <- rbind(biomassa_bacia, c(Bacias[i], comunidade))
}

rownames(biomassa_bacia) <- biomassa_bacia[,1]
biomassa_bacia <- as.data.frame(biomassa_bacia[,-1])

biomassa_bacia <- apply(biomassa_bacia, 2, as.numeric); rownames(biomassa_bacia) <- Bacias
biomassa_bacia <- as.data.frame(biomassa_bacia)

com_biomass_sum <- biomassa_bacia

########################################


write.csv(com, "data/com_por_trecho.csv")
write.csv(com_sum, "data/com_por_bacia.csv")
write.csv(com_biomass, "data/com_biomassa_por_trecho.csv")
write.csv(com_biomass_sum, "data/com_biomassa_por_bacia.csv")

