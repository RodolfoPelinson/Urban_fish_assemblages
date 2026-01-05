###Usando mv_abund

install.packages("mvabund")
library(mvabund)

#matriz resposta
assembleia_peixes_rm <- remove_sp(assembleia_peixes, 2)


#grupos de preditoras devem estar organizadas em listas (com nome)
predictors <- list(estrutura = data.frame(estrutura_PCs[,1:5]),
                   bacia = data.frame(bacia_PCs[,1:3]),
                   agua = data.frame(agua_PCs[,1:2]),
                   urbanizacao = decostand(data.frame(urb), "stand"))

#Usando a função
#Argumento resp é a matriz resposta, deve ser abundância! Sem transformar!!!!
#Pode usar mais argumentos que são passados para a função manyglm do pacote mvabund, ex: family (por padrão usa binomial negativa)
varpart_peixes <-varpart_manyglm(resp = assembleia_peixes_rm, pred = predictors, DF_adj_r2 = FALSE)

#Alguns resultados mais interessantes
varpart_peixes$R2_fractions_com
varpart_peixes$R2_fractions_sp$R2_full_fraction
round(varpart_peixes$R2_fractions_sp$R2_pure_fraction, 5)

#Teste para estrutura
anova_estrutura <- anova(varpart_peixes$model_null, varpart_peixes$models$estrutura, test = "LR", nBoot = 1000)
anova_estrutura

#Teste para estrutura pura
anova_estrutura_pura <- anova(varpart_peixes$models$`bacia-agua-urbanizacao`, varpart_peixes$models$`estrutura-bacia-agua-urbanizacao`, test = "LR", nBoot = 1000)
anova_estrutura_pura

#Olhando para os coeficientes do modelo só com estrutura
coefplot.manyglm(varpart_peixes$models$estrutura)

#Olhando para os coeficientes do modelo completo
coefplot.manyglm(varpart_peixes$models$`estrutura-bacia-agua-urbanizacao`)



######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################

devtools::install_github("JenniNiku/gllvm")
library(gllvm)
library(vegan)

#Preditores precisam estar todos em um mesmo data frame
predictors_gllvm <- data.frame(data.frame(estrutura_PCs[,1:2]),
                         data.frame(bacia_PCs[,1:2]),
                         data.frame(agua_PCs[,1:2]),
                         urb)

predictors_gllvm <- decostand(predictors_gllvm, method = "stand")

#Precisa usar fórmula, fiz isso aqui só pra automatizar caso queiramos mudar os preditores
formula_full <- formula( paste( "~", paste(colnames(predictors_gllvm), collapse = " + ") ) )

#Isso aqui efetivamente ajusta o modelo
model_gllvm_full <- gllvm(y = assembleia_peixes_rm, X = predictors_gllvm,formula = formula_full,
                      family = "negative.binomial",
                      num.lv = 1, method = "LA", control.start = list(n.init = 10, starting.val = "res"))

#Isso aqui é importante para fazer a partição de variancia, um vetor com os nomes das variaveis, e um com os numeros referentes ao grupo que fazem parte.
#Deve seguir a ordem com que as variaveis aparecem na matriz de preditores
groups <-c(rep(1, 2), rep(2, 2),rep(3, 2), rep(4, 1))
groups_names <-c("Estrutura", "Bacia","Água", "Urbanizacao")

#Fazer a partição de variância de fato
VP <- varPartitioning(
  model_gllvm_full,
  group = groups,
  groupnames = groups_names,
  adj.cov = TRUE,
  grouplvs = FALSE)

#Plotar a partição de variância.
plotVP(VP, col = c("violetred3", "#CDCD00", "green4", "azure4", "paleturquoise1"),
       las = 2, xlab = "", args.legend = list(), mar = c(10,4,2,.1))


