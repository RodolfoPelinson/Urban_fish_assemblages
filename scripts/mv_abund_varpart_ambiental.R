dir<-("C:/Users/rodol/OneDrive/repos/Urban_fish_assemblages")

source(paste(sep = "/",dir,"functions/remove_sp.R"))
source(paste(sep = "/",dir,"functions/R2_manyglm.R"))
source(paste(sep = "/",dir,"functions/forward_sel_manyglm.R"))

library(mvabund)
library(vegan)


assembleia_peixes <- read.csv(paste(sep = "/",dir,"data/com_por_bacia.csv"), row.names = 1)
assembleia_peixes <- assembleia_peixes[,-c(4,11)]
assembleia_peixes_rm <- remove_sp(assembleia_peixes, 1)

singletons_doubletons <- remove_sp(assembleia_peixes, 2, less_equal = TRUE)
singletons <- remove_sp(assembleia_peixes, 1, less_equal = TRUE)

sing_doub_ab <- rowSums(singletons_doubletons)
sing_ab <- rowSums(singletons)

sing_doub <- rowSums(decostand(singletons_doubletons, method = "pa")) 
sing <- rowSums(decostand(singletons, method = "pa")) 


assembleia_peixes_rm <- data.frame(assembleia_peixes_rm, Singletons  = sing_ab)


agua_PCs <- read.csv(paste(sep = "/",dir,"data/pcas_amb/agua_PCs.csv"), row.names = 1)
estrutura_PCs <- read.csv(paste(sep = "/",dir,"data/pcas_amb/estrutura_PCs.csv"), row.names = 1)
bacia_PCs <- read.csv(paste(sep = "/",dir,"data/pcas_amb/bacia_PCs.csv"), row.names = 1)

delineamento <- read.csv(paste(sep = "/",dir,"data/delineamento.csv"))

ncol(agua_PCs)
ncol(estrutura_PCs)
ncol(bacia_PCs)


############################################# FORWARD SELECTION ##########################################

est_FS <- forward_sel_manyglm(y = assembleia_peixes_rm, x = data.frame(est_PC1 = estrutura_PCs[,1],
                                                                       est_PC2 = estrutura_PCs[,2]), nBoot=999, quad = TRUE) #considerarei apenas os 5 primeiros eixos

agua_FS <- forward_sel_manyglm(y = assembleia_peixes_rm, x = data.frame(agua_PC1 = agua_PCs[,1],
                                                                        agua_PC2 = agua_PCs[,2]), nBoot=999, quad = TRUE) #considerarei apenas os 5 primeiros eixos

bacia_FS <- forward_sel_manyglm(y = assembleia_peixes_rm, x = data.frame(bacia_PC1 = bacia_PCs[,1],
                                                                         bacia_PC2 = bacia_PCs[,2]), nBoot=999, quad = TRUE) #considerarei apenas os 5 primeiros eixos


est_FS$new_x
agua_FS$new_x
bacia_FS$new_x

######################################################################################################################
######################################################################################################################
######################################## VARPART WITH FORWARD SELECTION ##############################################

predictors <- list(estrutura = est_FS$new_x,
                   agua = agua_FS$new_x,
                   bacia = bacia_FS$new_x)


varpart_peixes <- varpart_manyglm(resp = assembleia_peixes_rm, pred = predictors, DF_adj_r2 = TRUE)
varpart_peixes$R2_fractions_com
round(varpart_peixes$R2_fractions_sp$R2_full_fraction,4)
round(varpart_peixes$R2_fractions_sp$R2_pure_fraction,4)


full_model_sp <- varpart_peixes$R2_models_sp$estrutura.agua.bacia
names(full_model_sp) <- rownames(varpart_peixes$R2_models_sp)
ord_sp <- order(full_model_sp, decreasing = TRUE)

full_model_sp<- full_model_sp[ord_sp]



######PLOTS species
full_est <- varpart_peixes$R2_fractions_sp$R2_full_fraction$estrutura
names(full_est) <- rownames(varpart_peixes$R2_fractions_sp$R2_full_fraction)

full_agua <- varpart_peixes$R2_fractions_sp$R2_full_fraction$agua
names(full_agua) <- rownames(varpart_peixes$R2_fractions_sp$R2_full_fraction)

full_bacia <- varpart_peixes$R2_fractions_sp$R2_full_fraction$bacia
names(full_bacia) <- rownames(varpart_peixes$R2_fractions_sp$R2_full_fraction)


pure_est <- varpart_peixes$R2_fractions_sp$R2_pure_fraction$estrutura
names(pure_est) <- rownames(varpart_peixes$R2_fractions_sp$R2_pure_fraction)

pure_agua <- varpart_peixes$R2_fractions_sp$R2_pure_fraction$agua
names(pure_agua) <- rownames(varpart_peixes$R2_fractions_sp$R2_pure_fraction)

pure_bacia <- varpart_peixes$R2_fractions_sp$R2_pure_fraction$bacia
names(pure_bacia) <- rownames(varpart_peixes$R2_fractions_sp$R2_pure_fraction)


full_est<- full_est[ord_sp]
full_agua<- full_agua[ord_sp]
full_bacia<- full_bacia[ord_sp]

pure_est<- pure_est[ord_sp]
pure_agua<- pure_agua[ord_sp]
pure_bacia<- pure_bacia[ord_sp]

pure_est[pure_est < 0] <- 0
pure_agua[pure_agua < 0] <- 0
pure_bacia[pure_bacia < 0] <- 0


pure_frac_summed <- pure_est + pure_agua + pure_bacia

scale <- full_model_sp / pure_frac_summed

scale[scale > 1] <- 1

pure_frac_summed_scaled <- pure_frac_summed * scale

cbind(pure_frac_summed_scaled, full_model_sp)

pure_est_scaled  <- pure_est * scale
pure_agua_scaled  <- pure_agua * scale
pure_bacia_scaled  <- pure_bacia * scale



#########################3

full_model <- varpart_peixes$R2_models$`estrutura-agua-bacia`

pure_summed <- sum(varpart_peixes$R2_fractions_com$R2_pure_fraction)

scale_com <- full_model / pure_summed

pure_comm_scaled <- as.matrix(varpart_peixes$R2_fractions_com$R2_pure_fraction * scale_com)


##########################3


library(yarrr)


pdf(file = "C:/Users/rodol/OneDrive/repos/Urban_fish_assemblages/plots/varpart.pdf", width = 6, height = 5, pointsize = 12)
close.screen(all.screens = TRUE)
split.screen(matrix(c(0,0.3,0,1,
                      0.3,1,0,1), ncol = 4, nrow = 2, byrow = TRUE))

screen(2)
par(mar = c(12,1,1,1))
barplot(full_model_sp, ylim = c(0,1), las = 2, border = "white", col = "#C7C7C7", xaxt = "n", xaxs = "i", width = 1, space = 0.5, xlim = c(0,20), yaxt = "n")
par(new = TRUE)
pure_scaleds <- rbind(pure_est_scaled, pure_agua_scaled, pure_bacia_scaled)
barplot(pure_scaleds, ylim = c(0,1), las = 2, col = c("#98DF8A", "#9EDAE5", "#FFBB78"), border = "white", xaxt = "n", xaxs = "i", width = 1, space = 0.5, xlim = c(0,20), yaxt = "n")
box(bty = "l")
names <- colnames(pure_scaleds)
names <- gsub("_"," ", names)
axis(1, at = at_generator(first = 1, spacing = 1.5, n = length(names)), gap.axis = -10, tick = FALSE, labels = names, las = 2, font = 3, line = -0.5)
#title(ylab = "Likelihood ratio R²", cex.lab = 1.25)
axis(2, las = 2, line = 0, labels = FALSE)
par(new = TRUE, mar = c(0,0,0,0), bty = "n")
plot(NA, xlim = c(0,100), ylim = c(0,100), xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", )

legend(x = 99, y = 99, xjust = 1, yjust = 1, fill = c("#98DF8A", "#9EDAE5", "#FFBB78","#C7C7C7"),
       legend = c("Stream structure", "Water parameters", "Watershed", "Shared"), border = "white", bty = "n")


screen(1)
par(mar = c(12,4,1,1))
barplot(full_model, ylim = c(0,1), las = 2, border = "white", col = "#C7C7C7", xaxt = "n", xaxs = "i", width = 1, space = 0.5, xlim = c(0,2), yaxt = "n")
par(new = TRUE)

barplot(pure_comm_scaled, ylim = c(0,1), las = 2, col = c("#98DF8A", "#9EDAE5", "#FFBB78"), border = "white", xaxt = "n", xaxs = "i", width = 1, space = 0.5, xlim = c(0,2), yaxt = "n")
box(bty = "l")
axis(1, at = c(1), gap.axis = -10, tick = FALSE, labels = "Community", las = 1, font = 1, line = -0.5)
#title(ylab = "Likelihood ratio R²", cex.lab = 1.25)

title(ylab = "Likelihood ratio R²", cex.lab = 1.25)
axis(2, las = 2, line = 0)



dev.off()






pdf(file = "C:/Users/rodol/OneDrive/repos/Urban_fish_assemblages/plots/varpart_full.pdf", width = 6, height = 4.5, pointsize = 10)
close.screen(all.screens = TRUE)
split.screen(matrix(c(0,0.3,0,1,
                      0.3,1,0,1), ncol = 4, nrow = 2, byrow = TRUE))

screen(2)
par(mar = c(12,1,1,1))
full_sp <- rbind(full_est, full_agua, full_bacia)

barplot(full_sp, ylim = c(0,1), las = 2, col = c("#98DF8A", "#9EDAE5", "#FFBB78"),
        border = "transparent", xaxt = "n", xaxs = "i", width = 1, space = c(0,2), xlim = c(0,67), yaxt = "n", beside = TRUE)

box(bty = "l")
names <- colnames(pure_scaleds)
names <- gsub("_"," ", names)
#axis(1, at = at_generator(first = 3.5, spacing = 5, n = 12), gap.axis = -10, tick = T, labels = names, las = 2, font = 3, line = 0)
axis(1, at = at_generator(first = 3.5, spacing = 5, n = length(names)), gap.axis = -10, tick = T, labels = names, las = 2, font = 3, line = 0)
#title(ylab = "Likelihood ratio R²", cex.lab = 1.25)
axis(2, las = 2, line = 0, labels = FALSE)
par(new = TRUE, mar = c(0,0,0,0), bty = "n")
plot(NA, xlim = c(0,100), ylim = c(0,100), xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", )

legend(x = 99, y = 99, xjust = 1, yjust = 1, fill = c("#98DF8A", "#9EDAE5", "#FFBB78"),
       legend = c("Stream structure", "Water parameters", "Watershed"), border = "transparent", bty = "n")



screen(1)
par(mar = c(12,4,1,1))

full_comm <- as.matrix(varpart_peixes$R2_fractions_com$R2_full_fraction)


barplot(full_comm, ylim = c(0,1), las = 2, col = c("#98DF8A", "#9EDAE5", "#FFBB78"),
        border = "transparent", xaxt = "n", xaxs = "i", width = 1, space = c(0,1), xlim = c(0,5), yaxt = "n", beside = TRUE)
box(bty = "l")
axis(1, at = c(2.5), gap.axis = -10, tick = FALSE, labels = "Community", las = 1, font = 1, line = -0.5)
#title(ylab = "Likelihood ratio R²", cex.lab = 1.25)

title(ylab = "Likelihood ratio R²", cex.lab = 1.25)
axis(2, las = 2, line = 0)



dev.off()







names(full_model_sp) <- rownames(varpart_peixes$R2_fractions_sp$R2_full_fraction)


pure_est <- varpart_peixes$R2_models$`estrutura-agua-bacia` - varpart_peixes$R2_models$`agua-bacia`
pure_bacia <- varpart_peixes$R2_models$`estrutura-agua-bacia` - varpart_peixes$R2_models$`estrutura-agua`
pure_agua <- varpart_peixes$R2_models$`estrutura-agua-bacia` - varpart_peixes$R2_models$`estrutura-bacia`

#shared_baciaan_est <- varpart_peixes$R2_models$`estrutura-bacia` - (pure_est + pure_bacia)
#shared_baciaan_agua <- varpart_peixes$R2_models$`bacia-agua` - (pure_agua + pure_bacia)
#shared_est_agua <- varpart_peixes$R2_models$`estrutura-agua` - (pure_agua + pure_est)

shared_est_agua_sem_bacia <- full_model - (varpart_peixes$R2_models$bacia + pure_est + pure_agua)
shared_est_bacia_sem_agua <- full_model - (varpart_peixes$R2_models$agua + pure_est + pure_bacia)
shared_agua_bacia_sem_est <- full_model - (varpart_peixes$R2_models$estrutura + pure_agua + pure_bacia)

shared_all <- full_model - (pure_est + pure_bacia + pure_agua + shared_est_agua_sem_bacia + shared_est_bacia_sem_agua + shared_agua_bacia_sem_est)


unexplained<- 1 - full_model

p_est <- anova(varpart_peixes$model_null,
               varpart_peixes$models$estrutura,nBoot=999)
p_bacia <- anova(varpart_peixes$model_null,
                 varpart_peixes$models$bacia,nBoot=999)
p_agua <- anova(varpart_peixes$model_null,
                varpart_peixes$models$agua,nBoot=999)

p_est_pure <- anova(varpart_peixes$models$`agua-bacia`,
                    varpart_peixes$models$`estrutura-agua-bacia`,nBoot=999)
p_bacia_pure <- anova(varpart_peixes$models$`estrutura-agua`,
                      varpart_peixes$models$`estrutura-agua-bacia`,nBoot=999)
p_agua_pure <- anova(varpart_peixes$models$`estrutura-bacia`,
                     varpart_peixes$models$`estrutura-agua-bacia`,nBoot=999)

full_model <- anova(varpart_peixes$model_null,
                    varpart_peixes$models$`estrutura-agua-bacia`,nBoot=999)



#Predictions estanization
est_coefs <- varpart_peixes$models$estrutura$coefficients[2,]
est_IC_coefs <- varpart_peixes$models$estrutura$stderr.coefficients[2,] * qnorm(0.975)
est_upper_coefs <- est_coefs + est_IC_coefs
est_lower_coefs <- est_coefs - est_IC_coefs

names <- names(est_coefs)
names <- gsub("_"," ", names)

dev.off()
par(mar = c(4,12,1,1))
My_coefplot(mles = est_coefs, upper = est_upper_coefs,
            lower = est_lower_coefs, col_sig = "#98DF8A",
            cex_sig = 2, species_labels = names, yaxis_font = 3)
title(xlab = "Stream structure PC2", cex.lab = 1.5)

#newdata <- data.frame(estanizacao = seq(from = min(est), to = max(est), length.out = 100))
#predict.manyglm(varpart_peixes$models$estanizacao, newdata = newdata, type = "response")


#Predictions est_2anization
est_2_coefs <- varpart_peixes$models$estrutura$coefficients[3,]
est_2_IC_coefs <- varpart_peixes$models$estrutura$stderr.coefficients[3,] * qnorm(0.975)
est_2_upper_coefs <- est_2_coefs + est_2_IC_coefs
est_2_lower_coefs <- est_2_coefs - est_2_IC_coefs

names <- names(est_2_coefs)
names <- gsub("_"," ", names)

dev.off()
par(mar = c(4,12,1,1))
My_coefplot(mles = est_2_coefs, upper = est_2_upper_coefs,
            lower = est_2_lower_coefs, col_sig = "#98DF8A",
            cex_sig = 2, species_labels = names, yaxis_font = 3)
title(xlab = "(Stream structure PC2)²", cex.lab = 1.5)

#newdata <- data.frame(estanizacao = seq(from = min(est), to = max(est), length.out = 100))
#predict.manyglm(varpart_peixes$models$estanizacao, newdata = newdata, type = "response")




####################################################################################
### TRAIT



functional <- read.csv("data/functional_data.csv", row.names = 1)

functional[is.na(functional)] <- 0

rownames(functional)<- gsub(" ", "_", rownames(functional))


####Adding traits to singletons

#removendo Poecilia vivipara, que nao tem atributos funcionais
total_ab_singletons <- colSums(singletons)

functional_singletons <- functional[match(colnames(singletons_FN), rownames(functional)),]

functional_singleton <- colSums(functional_singletons * total_ab_singletons) / sum(total_ab_singletons)

functional <- rbind(functional, Singletons = functional_singleton)
###############################

assembleia_peixes_rm_funcional <- assembleia_peixes_rm[colnames(assembleia_peixes_rm) != "Poecilia_vivipara"]

functional <- functional[match(colnames(assembleia_peixes_rm_funcional), rownames(functional)),]



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



Model_pred_null <- traitglm(L = assembleia_peixes_rm_funcional, R = est_FS$new_x, Q = NULL)
Model_pred_trait_PC1 <- traitglm(L = assembleia_peixes_rm_funcional, R = est_FS$new_x, Q = data.frame(funcional_PCs[,1]))
Model_pred_trait_PC2 <- traitglm(L = assembleia_peixes_rm_funcional, R = est_FS$new_x, Q = data.frame(funcional_PCs[,2]))
Model_pred_trait_PC3 <- traitglm(L = assembleia_peixes_rm_funcional, R = est_FS$new_x, Q = data.frame(funcional_PCs[,3]))

anova_trait_PC1 <- anova(Model_pred_null, Model_pred_trait_PC1, nBoot = 9999)
anova_trait_PC2 <- anova(Model_pred_null, Model_pred_trait_PC2, nBoot = 9999)
anova_trait_PC3 <- anova(Model_pred_null, Model_pred_trait_PC3, nBoot = 9999)

c(Model_pred_trait$fourth.corner)
last <- nrow(Model_pred_trait$stderr.coefficients)
first <- (last - 2*1 ) + 1

Model_pred_trait$stderr.coefficients[first:last,]*qnorm(0.975)

#################################################################################

