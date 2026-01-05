#dieta <- read.csv("data/dieta_yulia.csv", sep = ";", dec = ",")
#write.csv(dieta, "data/dieta_yulia2.csv")
dieta <- read.csv("data/dieta_yulia2.csv")

index_dieta <- dieta[1,-c(1:2)]

dieta <- dieta[-1,]

sp<- dieta$Fish.ID
sp <- gsub("(", replacement = "", sp, fixed = TRUE)
sp <- gsub(")", replacement = "", sp, fixed = TRUE)
sp <- strsplit(sp, split = " ", fixed = TRUE)

especie <- rep(NA, length(sp))
individuo <- rep(NA, length(sp))
for(i in 1:length(sp)){
  especie[i] <- sp[[i]][2]
  individuo[i] <- sp[[i]][1]
}

dieta <- dieta[,-c(1,2)]

dieta <- apply(dieta, 2, as.numeric)
dieta <- apply(dieta, 2, round)
dieta_sum <- aggregate(dieta, by = list(as.factor(especie)), FUN = sum)
rownames(dieta_sum) <- dieta_sum$Group.1
dieta_sum <- dieta_sum[,-1]

diet_sum_index <- aggregate(t(dieta_sum), by = list(as.factor(unlist(c(index_dieta)))), FUN = sum)
rownames(diet_sum_index) <- diet_sum_index$Group.1
diet_sum_index <- t(diet_sum_index[,-1])





library(iNEXT)

inext_diet <- t(diet_sum_index)

inext <- iNEXT(x = inext_diet, q = c(0,1,2,3,4,5), datatype="abundance")
inext$iNextEst$coverage_based
inext$iNextEst$coverage_based[inext$iNextEst$coverage_based$Order.q == 2 & inext$iNextEst$coverage_based$SC > 0.999,]

inext <- iNEXT(x = inext_diet, q = c(0,1,2,3,4,5), datatype="abundance", size = 1300)

inext$iNextEst$size_based[inext$iNextEst$size_based$Order.q == 0 & inext$iNextEst$size_based$m == 1300,]
inext$iNextEst$size_based[inext$iNextEst$size_based$Order.q == 1 & inext$iNextEst$size_based$m == 1300,]
inext$iNextEst$size_based[inext$iNextEst$size_based$Order.q == 2 & inext$iNextEst$size_based$m == 1300,]
inext$iNextEst$size_based[inext$iNextEst$size_based$Order.q == 3 & inext$iNextEst$size_based$m == 1300,]

niche_amplitude <- inext$iNextEst$size_based[inext$iNextEst$size_based$Order.q == 0 & inext$iNextEst$size_based$m == 1300,][,c(1,3,8)]
niche_amplitude$q0 <- inext$iNextEst$size_based[inext$iNextEst$size_based$Order.q == 0 & inext$iNextEst$size_based$m == 1300,]$qD
niche_amplitude$q1 <- inext$iNextEst$size_based[inext$iNextEst$size_based$Order.q == 1 & inext$iNextEst$size_based$m == 1300,]$qD
niche_amplitude$q2 <- inext$iNextEst$size_based[inext$iNextEst$size_based$Order.q == 2 & inext$iNextEst$size_based$m == 1300,]$qD
niche_amplitude$q3 <- inext$iNextEst$size_based[inext$iNextEst$size_based$Order.q == 3 & inext$iNextEst$size_based$m == 1300,]$qD
niche_amplitude$q4 <- inext$iNextEst$size_based[inext$iNextEst$size_based$Order.q == 4 & inext$iNextEst$size_based$m == 1300,]$qD

write.csv(niche_amplitude, "data/niche_amplitude.csv")

