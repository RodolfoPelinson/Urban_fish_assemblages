dir<-("C:/Users/rodol/OneDrive/repos/Urban_fish_assemblages")

source(paste(sep = "/",dir,"functions/remove_sp.R"))
source(paste(sep = "/",dir,"functions/R2_manyglm.R"))
source(paste(sep = "/",dir,"functions/forward_sel_manyglm.R"))

library(mvabund)
library(vegan)


assembleia_peixes <- read.csv(paste(sep = "/",dir,"data/com_por_bacia.csv"), row.names = 1)
assembleia_peixes <- assembleia_peixes[,-c(4,11)]
assembleia_peixes_rm <- remove_sp(com = assembleia_peixes, n_sp = 1)

singletons_doubletons <- remove_sp(assembleia_peixes, 2, less_equal = TRUE)
singletons <- remove_sp(assembleia_peixes, 1, less_equal = TRUE)

sing_doub_ab <- rowSums(singletons_doubletons)
sing_ab <- rowSums(singletons)

assembleia_peixes_rm <- data.frame(assembleia_peixes_rm, Singletons = sing_ab)

agua_PCs <- read.csv(paste(sep = "/",dir,"data/pcas_amb/agua_PCs.csv"), row.names = 1)
estrutura_PCs <- read.csv(paste(sep = "/",dir,"data/pcas_amb/estrutura_PCs.csv"), row.names = 1)
bacia_PCs <- read.csv(paste(sep = "/",dir,"data/pcas_amb/bacia_PCs.csv"), row.names = 1)

delineamento <- read.csv(paste(sep = "/",dir,"data/delineamento.csv"))

urb <- data.frame(urb = delineamento$urbana)

urb <- decostand(urb, method = "stand")
agua_PCs <- decostand(agua_PCs, method = "stand")
estrutura_PCs <- decostand(estrutura_PCs, method = "stand")
bacia_PCs <- decostand(bacia_PCs, method = "stand")


ncol(agua_PCs)
ncol(estrutura_PCs)
ncol(bacia_PCs)


library(ade4)
library(adespatial)

dist_euclid <- read.csv(paste(sep = "/",dir,"data/dist/Matriz_distancia_matriz_euclidiana.csv"), row.names = 1)
dist_euclid <- as.dist(dist_euclid)

dbmem_euclid <- dbmem(dist_euclid, thresh = NULL, MEM.autocor = c("positive", "non-null", "all", "negative"), store.listw = TRUE, silent = FALSE)

dbmem_euclid <- decostand(dbmem_euclid, method = "stand")

############################################# FORWARD SELECTION ##########################################

env_data.frame <- data.frame(est_PC1 = estrutura_PCs[,1],
                             est_PC2 = estrutura_PCs[,2],
                             est_PC3 = estrutura_PCs[,3],
                             est_PC4 = estrutura_PCs[,4],
                             est_PC5 = estrutura_PCs[,5],
                             agua_PC1 = agua_PCs[,1],
                             agua_PC2 = agua_PCs[,2],
                             agua_PC3 = agua_PCs[,3],
                             agua_PC4 = agua_PCs[,4],
                             agua_PC5 = agua_PCs[,5],
                             bacia_PC1 = bacia_PCs[,1],
                             bacia_PC2 = bacia_PCs[,2],
                             bacia_PC3 = bacia_PCs[,3],
                             bacia_PC4 = bacia_PCs[,4],
                             bacia_PC5 = bacia_PCs[,5]) #considerarei apenas os 5 primeiros eixos

library(corrplot)
corrplot(cor(env_data.frame))

amb_FS <- forward_sel_manyglm(y = assembleia_peixes_rm, x = env_data.frame, nBoot=999, quad = TRUE) #considerarei apenas os 5 primeiros eixos

esp_FS <- forward_sel_manyglm(y = assembleia_peixes_rm, x = data.frame(dbmem_euclid), nBoot=999, quad = FALSE) 


assembleia_peixes_rm_mv <- mvabund(assembleia_peixes_rm)

model_null <- manyglm(assembleia_peixes_rm_mv ~ 1,data = urb)
model_urb <- manyglm(assembleia_peixes_rm_mv ~ urb,data = urb)
model_urb_quad <- manyglm(assembleia_peixes_rm_mv ~ urb + I(urb^2),data = urb)

anova(model_null, model_urb, model_urb_quad)
anova(model_null, model_urb_quad)

######################################################################################################################
######################################################################################################################
######################################## VARPART WITH FORWARD SELECTION ##############################################

predictors <- list(ambiente = amb_FS$new_x,
                   urbanizacao = data.frame(urbanizacao = urb$urb),
                   MEMs = esp_FS$new_x)

varpart_peixes <- varpart_manyglm(resp = assembleia_peixes_rm, pred = predictors, DF_adj_r2 = TRUE)

varpart_peixes$R2_fractions_com
round(varpart_peixes$R2_fractions_sp$R2_full_fraction,4)
round(varpart_peixes$R2_fractions_sp$R2_pure_fraction,4)


full_model_sp <- varpart_peixes$R2_models_sp$ambiente.urbanizacao.MEMs
names(full_model_sp) <- rownames(varpart_peixes$R2_models_sp)
ord_sp <- order(full_model_sp, decreasing = TRUE)

full_model_sp <- full_model_sp[ord_sp]







######PLOTS species
full_amb <- varpart_peixes$R2_fractions_sp$R2_full_fraction$ambiente
names(full_amb) <- rownames(varpart_peixes$R2_fractions_sp$R2_full_fraction)

full_urb <- varpart_peixes$R2_fractions_sp$R2_full_fraction$urbanizacao
names(full_urb) <- rownames(varpart_peixes$R2_fractions_sp$R2_full_fraction)

full_MEMs <- varpart_peixes$R2_fractions_sp$R2_full_fraction$MEMs
names(full_MEMs) <- rownames(varpart_peixes$R2_fractions_sp$R2_full_fraction)


pure_amb <- varpart_peixes$R2_fractions_sp$R2_pure_fraction$ambiente
names(pure_amb) <- rownames(varpart_peixes$R2_fractions_sp$R2_pure_fraction)

pure_urb <- varpart_peixes$R2_fractions_sp$R2_pure_fraction$urbanizacao
names(pure_urb) <- rownames(varpart_peixes$R2_fractions_sp$R2_pure_fraction)

pure_MEMs <- varpart_peixes$R2_fractions_sp$R2_pure_fraction$MEMs
names(pure_MEMs) <- rownames(varpart_peixes$R2_fractions_sp$R2_pure_fraction)


full_amb<- full_amb[ord_sp]
full_urb<- full_urb[ord_sp]
full_MEMs<- full_MEMs[ord_sp]

pure_amb<- pure_amb[ord_sp]
pure_urb<- pure_urb[ord_sp]
pure_MEMs<- pure_MEMs[ord_sp]

pure_amb[pure_amb < 0] <- 0
pure_urb[pure_urb < 0] <- 0
pure_MEMs[pure_MEMs < 0] <- 0


pure_frac_summed <- pure_amb + pure_urb + pure_MEMs

scale <- full_model_sp / pure_frac_summed

scale[scale > 1] <- 1

pure_frac_summed_scaled <- pure_frac_summed * scale

cbind(pure_frac_summed_scaled, full_model_sp)

pure_amb_scaled  <- pure_amb * scale
pure_urb_scaled  <- pure_urb * scale
pure_MEMs_scaled  <- pure_MEMs * scale




#########################3

varpart_peixes$R2_fractions_com

full_model <- varpart_peixes$R2_models$`ambiente-urbanizacao-MEMs`

pure_summed <- sum(varpart_peixes$R2_fractions_com$R2_pure_fraction)

pure_comm <- as.matrix(varpart_peixes$R2_fractions_com$R2_pure_fraction)

pure_comm[pure_comm < 0] <- 0

full_com <- as.matrix(varpart_peixes$R2_fractions_com$R2_full_fraction)

##########################3




library(yarrr)


pdf(file = "C:/Users/rodol/OneDrive/repos/Urban_fish_assemblages/plots/varpart2.pdf", width = 6, height = 4.5, pointsize = 10)
close.screen(all.screens = TRUE)
split.screen(matrix(c(0,0.3,0,1,
                      0.3,1,0,1), ncol = 4, nrow = 2, byrow = TRUE))

screen(2)
par(mar = c(12,1,1,1))
barplot(full_model_sp, ylim = c(0,1), las = 2, border = "transparent", col = "#C7C7C7", xaxt = "n", xaxs = "i", width = 1, space = 0.5, xlim = c(0,18.5), yaxt = "n")
par(new = TRUE)
pure_scaleds <- rbind(pure_amb_scaled, pure_urb_scaled, pure_MEMs_scaled)
barplot(pure_scaleds, ylim = c(0,1), las = 2, col = c("#DBDB8D", "#AEC7E8", "#C49C94"), border = "transparent", xaxt = "n", xaxs = "i", width = 1, space = 0.5, xlim = c(0,18.5), yaxt = "n")
box(bty = "l")
names <- colnames(pure_scaleds)
names <- gsub("_"," ", names)
axis(1, at = at_generator(first = 1, spacing = 1.5, n = 7), gap.axis = -10, tick = TRUE, labels = names, las = 2, font = 3, line = 0)
#title(ylab = "Likelihood ratio R²", cex.lab = 1.25)
axis(2, las = 2, line = 0, labels = FALSE)
par(new = TRUE, mar = c(0,0,0,0), bty = "n")
plot(NA, xlim = c(0,100), ylim = c(0,100), xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", )

legend(x = 99, y = 99, xjust = 1, yjust = 1, fill = c("#DBDB8D", "#AEC7E8", "#C49C94","#C7C7C7"),
       legend = c("Environmental parameters", "Urban cover", "Spatial filters", "Shared"), border = "transparent", bty = "n")



screen(1)
par(mar = c(12,4,1,1))
barplot(full_model, ylim = c(0,1), las = 2, border = "transparent", col = "#C7C7C7", xaxt = "n", xaxs = "i", width = 1, space = 0.5, xlim = c(0,2), yaxt = "n")
par(new = TRUE)

barplot(pure_comm, ylim = c(0,1), las = 2, col = c("#DBDB8D", "#AEC7E8", "#C49C94"), border = "transparent", xaxt = "n", xaxs = "i", width = 1, space = 0.5, xlim = c(0,2), yaxt = "n")
box(bty = "l")
axis(1, at = c(1), gap.axis = -10, tick = FALSE, labels = "Community", las = 1, font = 1, line = -0.5)
#title(ylab = "Likelihood ratio R²", cex.lab = 1.25)

title(ylab = "Likelihood ratio R²", cex.lab = 1.25)
axis(2, las = 2, line = 0)



dev.off()







pdf(file = "C:/Users/rodol/OneDrive/repos/Urban_fish_assemblages/plots/varpart2_full.pdf", width = 6, height = 4.5, pointsize = 10)
close.screen(all.screens = TRUE)
split.screen(matrix(c(0,0.3,0,1,
                      0.3,1,0,1), ncol = 4, nrow = 2, byrow = TRUE))

screen(2)
par(mar = c(12,1,1,1))
full_sp <- rbind(full_amb, full_urb, full_MEMs)

barplot(full_sp, ylim = c(0,1), las = 2, col = c("#DBDB8D", "#AEC7E8", "#C49C94"),
        border = "transparent", xaxt = "n", xaxs = "i", width = 1, space = c(0,2), xlim = c(0,62), yaxt = "n", beside = TRUE)

box(bty = "l")
names <- colnames(pure_scaleds)
names <- gsub("_"," ", names)
#axis(1, at = at_generator(first = 3.5, spacing = 5, n = 12), gap.axis = -10, tick = T, labels = names, las = 2, font = 3, line = 0)
axis(1, at = at_generator(first = 3.5, spacing = 5, n = 12), gap.axis = -10, tick = T, labels = names, las = 2, font = 3, line = 0)
#title(ylab = "Likelihood ratio R²", cex.lab = 1.25)
axis(2, las = 2, line = 0, labels = FALSE)
par(new = TRUE, mar = c(0,0,0,0), bty = "n")
plot(NA, xlim = c(0,100), ylim = c(0,100), xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", )

legend(x = 99, y = 99, xjust = 1, yjust = 1, fill = c("#DBDB8D", "#AEC7E8", "#C49C94"),
       legend = c("Environmental parameters", "Urban cover", "Spatial filters"), border = "transparent", bty = "n")



screen(1)
par(mar = c(12,4,1,1))

full_comm <- as.matrix(varpart_peixes$R2_fractions_com$R2_full_fraction)


barplot(full_comm, ylim = c(0,1), las = 2, col = c("#DBDB8D", "#AEC7E8", "#C49C94"),
        border = "transparent", xaxt = "n", xaxs = "i", width = 1, space = c(0,1), xlim = c(0,5), yaxt = "n", beside = TRUE)
box(bty = "l")
axis(1, at = c(2.5), gap.axis = -10, tick = FALSE, labels = "Community", las = 1, font = 1, line = -0.5)
#title(ylab = "Likelihood ratio R²", cex.lab = 1.25)

title(ylab = "Likelihood ratio R²", cex.lab = 1.25)
axis(2, las = 2, line = 0)



dev.off()














unexplained<- 1 - full_model

p_amb <- anova(varpart_peixes$model_null,
               varpart_peixes$models$ambiente,nBoot=999)
p_MEMs <- anova(varpart_peixes$model_null,
                 varpart_peixes$models$MEMs,nBoot=999)
p_urbanizacao <- anova(varpart_peixes$model_null,
                varpart_peixes$models$urbanizacao,nBoot=999)

varpart_peixes$R2_fractions_com

pchisq(31.64, 2)

p_amb_pure <- anova(varpart_peixes$models$`urbanizacao-MEMs`,
                    varpart_peixes$models$`ambiente-urbanizacao-MEMs`,nBoot=999)
p_MEMs_pure <- anova(varpart_peixes$models$`ambiente-urbanizacao`,
                      varpart_peixes$models$`ambiente-urbanizacao-MEMs`,nBoot=999)
p_urbanizacao_pure <- anova(varpart_peixes$models$`ambiente-MEMs`,
                     varpart_peixes$models$`ambiente-urbanizacao-MEMs`,nBoot=999)

full_model <- anova(varpart_peixes$model_null,
                    varpart_peixes$models$`ambiente-urbanizacao-MEMs`,nBoot=999)








####################################################################################
### TRAIT



functional <- read.csv("data/functional_data.csv", row.names = 1)

rownames(functional)<- gsub(" ", "_", rownames(functional))

functional <- functional[match(colnames(assembleia_peixes_rm), rownames(functional)),]

functional[is.na(functional)] <- 0

library(vegan)

cor(functional)


functional_stand <- decostand(functional, method = "stand")

pca_funcional <- rda(functional_stand)

importance_funcional <- round(pca_funcional$CA$eig/sum(pca_funcional$CA$eig),2)
Eigenvalues_funcional <- data.frame(autovalores = pca_funcional$CA$eig,
                                    importance = importance_funcional)



sum(importance_funcional[1:3])
funcional_PCs <- pca_funcional$CA$u
funcional_loadings <- pca_funcional$CA$v

funcional_loadings_filtrados <- funcional_loadings[which(funcional_loadings[,1] > 0.2 | funcional_loadings[,1] < -0.2 |
                                                           funcional_loadings[,2] > 0.2 | funcional_loadings[,2] < -0.2 | 
                                                           funcional_loadings[,3] > 0.2 | funcional_loadings[,3] < -0.2),1:3]



Model_pred_null <- traitglm(L = assembleia_peixes_rm, R = urb, Q = NULL)
Model_pred_trait <- traitglm(L = assembleia_peixes_rm, R = urb, Q = data.frame(funcional_PCs[,1:3]))

anova_trait <- anova(Model_pred_null, Model_pred_trait, nBoot = 999)

c(Model_pred_trait$fourth.corner)
Model_pred_trait$stderr.coefficients[10:12,]*qnorm(0.975)

#################################################################################

#Predictions


urb_coefs <- varpart_peixes$models$urbanizacao$coefficients[2,]
urb_IC_coefs <- varpart_peixes$models$urbanizacao$stderr.coefficients[2,] * qnorm(0.975)

coefplot.manyglm(varpart_peixes$models$urbanizacao)
coefplot.manyglm(varpart_peixes$models$ambiente)
#coefplot.manyglm(varpart_peixes$models$MEMs)


newdata <- data.frame(urbanizacao = seq(from = min(urb), to = max(urb), length.out = 100))
predict.manyglm(varpart_peixes$models$urbanizacao, newdata = newdata, type = "response")





######################################################################


library(glmmTMB)
library(DHARMa)
library(car)


###Abundances of singletons

sing_mod_null_ab <- glmmTMB(sing_ab ~ 1, family = "nbinom2", data = data.frame(amb_FS$new_x, urb = urb$urb))
sing_mod_urb_ab <- glmmTMB(sing_ab ~ urb, family = "nbinom2", data = data.frame(amb_FS$new_x, urb = urb$urb))
sing_mod_urb_quad_ab <- glmmTMB(sing_ab ~ urb + I(urb^2), family = "nbinom2", data = data.frame(amb_FS$new_x, urb = urb$urb))
sing_mod_amb_ab <- glmmTMB(sing_ab ~ est_PC2, family = "nbinom2", data = data.frame(amb_FS$new_x, urb = urb$urb))
sing_mod_amb_quad_ab <- glmmTMB(sing_ab ~ est_PC2 + I(est_PC2^2), family = "nbinom2", data = data.frame(amb_FS$new_x, urb = urb$urb))

sing_mod_urb_amb_ab <- glmmTMB(sing_ab ~ urb + est_PC2, family = "nbinom2", data = data.frame(amb_FS$new_x, urb = urb$urb))
sing_mod_urb_quad_amb_ab <- glmmTMB(sing_ab ~ urb + I(urb^2) + est_PC2, family = "nbinom2", data = data.frame(amb_FS$new_x, urb = urb$urb))

sing_mod_urb_amb_quad_ab <- glmmTMB(sing_ab ~ urb + est_PC2+ I(est_PC2^2), family = "nbinom2", data = data.frame(amb_FS$new_x, urb = urb$urb))
sing_mod_urb_quad_amb_quad_ab <- glmmTMB(sing_ab ~ urb + I(urb^2) + est_PC2+ I(est_PC2^2), family = "nbinom2", data = data.frame(amb_FS$new_x, urb = urb$urb))


#First, we test the hypothesis of a linear effect of urbanization
anova(sing_mod_null_ab, sing_mod_urb_ab) #NO

#Second, we test the hypothesis of a quadratic effect of urbanization
anova(sing_mod_null_ab, sing_mod_urb_quad_ab) #NO

#Now we test the hypothesis of a linear effect of environment
anova(sing_mod_null_ab, sing_mod_amb_ab) #YES!!!!!!!

#Now we test the hypothesis of quadratic effect of environment
anova(sing_mod_null_ab, sing_mod_amb_quad_ab) #Almost!!!!!!

#We also wanna test the same urbanization hypothesis after accounting for environmental effects.
anova(sing_mod_amb_ab, sing_mod_urb_amb_ab) #NO

#We also wanna test the same quadratic urbanization hypothesis after accounting for environmental effects.
anova(sing_mod_amb_ab, sing_mod_urb_quad_amb_ab) #NO

#We also wanna test the same urbanization hypothesis after accounting for quadratic environmental effects.
anova(sing_mod_amb_quad_ab, sing_mod_urb_amb_quad_ab) #NO

#We also wanna test the same quadratic urbanization hypothesis after accounting for quadratic environmental effects.
anova(sing_mod_amb_quad_ab, sing_mod_urb_quad_amb_quad_ab) #NO







###Abundances of singletons and doubletons

soub_mod_null_ab <- glmmTMB(sing_doub_ab ~ 1, family = "nbinom2", data = data.frame(amb_FS$new_x, urb = urb$urb))
soub_mod_urb_ab <- glmmTMB(sing_doub_ab ~ urb, family = "nbinom2", data = data.frame(amb_FS$new_x, urb = urb$urb))
soub_mod_urb_quad_ab <- glmmTMB(sing_doub_ab ~ urb + I(urb^2), family = "nbinom2", data = data.frame(amb_FS$new_x, urb = urb$urb))
soub_mod_amb_ab <- glmmTMB(sing_doub_ab ~ est_PC2, family = "nbinom2", data = data.frame(amb_FS$new_x, urb = urb$urb))
soub_mod_amb_quad_ab <- glmmTMB(sing_doub_ab ~ est_PC2 + I(est_PC2^2), family = "nbinom2", data = data.frame(amb_FS$new_x, urb = urb$urb))

soub_mod_urb_amb_ab <- glmmTMB(sing_doub_ab ~ urb + est_PC2, family = "nbinom2", data = data.frame(amb_FS$new_x, urb = urb$urb))
soub_mod_urb_quad_amb_ab <- glmmTMB(sing_doub_ab ~ urb + I(urb^2) + est_PC2, family = "nbinom2", data = data.frame(amb_FS$new_x, urb = urb$urb))

soub_mod_urb_amb_quad_ab <- glmmTMB(sing_doub_ab ~ urb + est_PC2+ I(est_PC2^2), family = "nbinom2", data = data.frame(amb_FS$new_x, urb = urb$urb))
soub_mod_urb_quad_amb_quad_ab <- glmmTMB(sing_doub_ab ~ urb + I(urb^2) + est_PC2+ I(est_PC2^2), family = "nbinom2", data = data.frame(amb_FS$new_x, urb = urb$urb))


#First, we test the hypothesis of a linear effect of urbanization
anova(soub_mod_null_ab, soub_mod_urb_ab) #NO

#Second, we test the hypothesis of a quadratic effect of urbanization
anova(soub_mod_null_ab, soub_mod_urb_quad_ab) #NO

#Now we test the hypothesis of a linear effect of environment
anova(soub_mod_null_ab, soub_mod_amb_ab) #NO

#Now we test the hypothesis of a quadratic effect of environment
anova(soub_mod_null_ab, soub_mod_amb_quad_ab) #YESS!!!!!!!!!!!

#We also wanna test the same urbanization hypothesis after accounting for environmental effects.
anova(soub_mod_amb_ab, soub_mod_urb_amb_ab) #NO

#We also wanna test the same quadratic urbanization hypothesis after accounting for environmental effects.
anova(soub_mod_amb_ab, soub_mod_urb_quad_amb_ab) #NO

#We also wanna test the same urbanization hypothesis after accounting for quadratic environmental effects.
anova(soub_mod_amb_quad_ab, soub_mod_urb_amb_quad_ab) #NO

#We also wanna test the same quadratic urbanization hypothesis after accounting for quadratic environmental effects.
anova(soub_mod_amb_quad_ab, soub_mod_urb_quad_amb_quad_ab) #NO




