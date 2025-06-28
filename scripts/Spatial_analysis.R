dist_euclid <- read.csv("data/dist/Matriz_distancia_matriz_euclidiana.csv", row.names = 1)
dist_euclid <- as.dist(dist_euclid)

dist_fluv <- read.csv("data/dist/Matriz_distancia_fluvial.csv", row.names = 1)
dist_fluv <- as.dist(dist_fluv)

dist_decl <- read.csv("data/dist/Matriz_distancia_declividade_v2.csv", row.names = 1)
dist_decl <- as.dist(dist_decl)


library(ade4)
library(adespatial)


dbmem_euclid <- dbmem(dist_euclid, thresh = NULL, MEM.autocor = c("positive", "non-null", "all", "negative"), store.listw = TRUE, silent = FALSE)

dbmem_fluv <- dbmem(dist_fluv, thresh = NULL, MEM.autocor = c("positive", "non-null", "all", "negative"), store.listw = TRUE, silent = FALSE)

dbmem_dec <- dbmem(dist_decl, thresh = NULL, MEM.autocor = c("positive", "non-null", "all", "negative"), store.listw = TRUE, silent = FALSE)

remove <- which(is.na(match(rownames(estrutura_PCs), rownames(dbmem_fluv))))

###Usando mv_abund

library(mvabund)
library(vegan)

#matriz resposta
assembleia_peixes_rm <- remove_sp(assembleia_peixes, 2)



estrutura_FS <- forward_sel_manyglm(y = assembleia_peixes_rm, x = data.frame(estrutura_PCs[,1:5]))
bacia_FS <- forward_sel_manyglm(y = assembleia_peixes_rm, x = data.frame(bacia_PCs[,1:3]))
agua_FS <- forward_sel_manyglm(y = assembleia_peixes_rm, x = data.frame(agua_PCs[,1:2]))
urbanizacao_FS <- forward_sel_manyglm(y = assembleia_peixes_rm, x = decostand(data.frame(urb), "stand"))

dbmem_euclid_FS <- forward_sel_manyglm(y = assembleia_peixes_rm, x = data.frame(dbmem_euclid))
dbmem_fluv_FS <- forward_sel_manyglm(y = assembleia_peixes_rm[-remove,], x = data.frame(dbmem_fluv))
dbmem_dec_FS <- forward_sel_manyglm(y = assembleia_peixes_rm, x = data.frame(dbmem_dec))


#grupos de preditoras devem estar organizadas em listas (com nome)
predictors_euclid <- list(estrutura = estrutura_FS$new_x,
                   bacia = bacia_FS$new_x,
                   agua = agua_FS$new_x,
                   #urbanizacao = urbanizacao_FS$new_x,
                   MEMs_euclid = dbmem_euclid_FS$new_x)

predictors_fluv <- list(estrutura = data.frame(estrutura_FS$new_x[-remove,]),
                          bacia = data.frame(bacia_FS$new_x[-remove,]),
                          agua = data.frame(agua_FS$new_x[-remove,]),
                          #urbanizacao = data.frame(urbanizacao_FS$new_x[-remove,]),
                          MEMs_fluv = dbmem_fluv_FS$new_x)

predictors_dec <- list(estrutura = estrutura_FS$new_x,
                          bacia = bacia_FS$new_x,
                          agua = agua_FS$new_x,
                          #urbanizacao = urbanizacao_FS$new_x,
                          MEMs_dec = dbmem_dec_FS$new_x)







#Usando a função
#Argumento resp é a matriz resposta, deve ser abundância! Sem transformar!!!!
#Pode usar mais argumentos que são passados para a função manyglm do pacote mvabund, ex: family (por padrão usa binomial negativa)
varpart_peixes_euclid <-varpart_manyglm(resp = assembleia_peixes_rm, pred = predictors_euclid, DF_adj_r2 = FALSE)
round(varpart_peixes_euclid$R2_fractions_com,7)
round(varpart_peixes_euclid$R2_fractions_sp$R2_full_fraction,7)
round(varpart_peixes_euclid$R2_fractions_sp$R2_pure_fraction,7)


varpart_peixes_fluv <-varpart_manyglm(resp = assembleia_peixes_rm[-remove,], pred = predictors_fluv, DF_adj_r2 = FALSE)
round(varpart_peixes_fluv$R2_fractions_com,7)
round(varpart_peixes_fluv$R2_fractions_sp$R2_full_fraction,7)
round(varpart_peixes_fluv$R2_fractions_sp$R2_pure_fraction,7)


varpart_peixes_dec <-varpart_manyglm(resp = assembleia_peixes_rm, pred = predictors_dec, DF_adj_r2 = FALSE)
round(varpart_peixes_dec$R2_fractions_com,7)
round(varpart_peixes_dec$R2_fractions_sp$R2_full_fraction,7)
round(varpart_peixes_dec$R2_fractions_sp$R2_pure_fraction,7)


##############################################################
##############################################################


assembleia_peixes_rm <- remove_sp(assembleia_peixes, 2)


estrutura_FSV <- forward_selection_vegan(decostand(assembleia_peixes_rm, method = "hel"), data.frame(estrutura_PCs[,1:5])) #PC1
bacia_FSV <-forward_selection_vegan(decostand(assembleia_peixes_rm, method = "hel"), data.frame(bacia_PCs[,1:3])) #PC1 e PC3
agua_FSV <-forward_selection_vegan(decostand(assembleia_peixes_rm, method = "hel"), data.frame(agua_PCs[,1:2])) #PC1
#urb_FSV <-forward_selection_vegan(decostand(assembleia_peixes_rm, method = "hel"), decostand(data.frame(urb), "stand")) #urb
dbmen_eucl_FSV <-forward_selection_vegan(decostand(assembleia_peixes_rm, method = "hel"), data.frame(dbmem_euclid)) #MEM5
dbmen_fluv_FSV <-forward_selection_vegan(decostand(assembleia_peixes_rm, method = "hel")[-remove,], data.frame(dbmem_fluv)) #Nenhuma
dbmen_dec_FSV <-forward_selection_vegan(decostand(assembleia_peixes_rm, method = "hel"), data.frame(dbmem_dec)) #MEM4


varpart_euclid <- varpart(assembleia_peixes_rm,
                          estrutura_FSV$selected_variables, #X1
                          bacia_FSV$selected_variables, #X2
                          agua_FSV$selected_variables, #X3
                          dbmen_eucl_FSV$selected_variables, #X4
                          transfo="hel")
plot(varpart_euclid, bg=2:4, Xnames = c("Estrutura", "Bacia", "Água", "MEM_euclid"))

anova(rda(assembleia_peixes_rm, dbmen_eucl_FSV$selected_variables))

anova(rda(assembleia_peixes_rm, dbmen_eucl_FSV$selected_variables, data.frame(bacia_FSV$selected_variables,
                                                                              agua_FSV$selected_variables,
                                                                              estrutura_FSV$selected_variables)))


varpart_fluv <- varpart(assembleia_peixes_rm[-remove,],
                          estrutura_FSV$selected_variables[-remove,], #X1
                          bacia_FSV$selected_variables[-remove,], #X2
                          agua_FSV$selected_variables[-remove,], #X3
                          dbmen_fluv_FSV$selected_variables, #X4
                          transfo="hel")
plot(varpart_fluv, bg=2:4, Xnames = c("Estrutura", "Bacia", "Água", "MEM_fluv"))



varpart_dec <- varpart(assembleia_peixes_rm,
                          estrutura_FSV$selected_variables, #X1
                          bacia_FSV$selected_variables, #X2
                          agua_FSV$selected_variables, #X3
                          dbmen_dec_FSV$selected_variables, #X4
                          transfo="hel")
plot(varpart_dec, bg=2:4, Xnames = c("Estrutura", "Bacia", "Água", "MEM_dec"))

anova(rda(assembleia_peixes_rm, dbmen_dec_FSV$selected_variables, data.frame(bacia_FSV$selected_variables,
                                                                              agua_FSV$selected_variables,
                                                                              estrutura_FSV$selected_variables)))

anova(rda(assembleia_peixes_rm, dbmen_dec_FSV$selected_variables))

      