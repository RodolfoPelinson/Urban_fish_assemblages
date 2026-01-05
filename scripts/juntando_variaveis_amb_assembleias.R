agua <- read.csv("data/dados_ambientais_agua_2021.csv")
bacia <- read.csv("data/dados_ambientais_bacia_3.csv")
bacia2 <- read.csv("data/dados_ambientais_bacia_2.csv")
estrutura <- read.csv("data/dados_ambientais_estrutura_por_riacho.csv")





##############################################################################



dados_peixes <- read.table("data/com_por_bacia.csv")

nrow(agua)
nrow(bacia)
nrow(estrutura)
nrow(dados_peixes)



#Removendo segunda linha da planilha de estrutura e adicionando a informação na primeiras
primeira <- colnames(estrutura)
segunda <- estrutura[1,]

nova_primeira <- paste(segunda, primeira, sep = "_")
nova_primeira <- gsub(x = nova_primeira, pattern = "NA_", replacement = "")
nova_primeira <- gsub(x = nova_primeira, pattern = ".1", replacement = "")
nova_primeira <- gsub(x = nova_primeira, pattern = ".2", replacement = "")
nova_primeira <- gsub(x = nova_primeira, pattern = ".3", replacement = "")

estrutura <- estrutura[-1,]
colnames(estrutura) <- nova_primeira

rownames(estrutura) <- estrutura$bacia

#Separando a variavel urbanização pra comparar com os valores da planilha bacia
urb <- estrutura$urbanizacao
names(urb) <- rownames(estrutura)

estrutura<- estrutura[,-c(1:2)]

nrow(estrutura)


#Removendo os coeficientes de variação
estrutura <- estrutura[colnames(estrutura) != "prof_cv" &
                         colnames(estrutura) != "larg_cv"&
                         colnames(estrutura) != "velocidade_cv"]



#Criando uma variavel alternativa para descarga
estrutura$descarga_vel_dim <- estrutura$velocidade_media * estrutura$larg_media * estrutura$prof_media 

############################################################################


##############################################################################
################ ESTIMATING MISSING VALUES IN THE WATER MATRIX ###############

#Removendo bacias não usadas e reordenando a planilha de água de acordo com a planilha de estrutura

posicoes <- match(rownames(estrutura), agua$Group.1)
na_posicoes <- which(is.na(posicoes))

new_bacias <- data.frame(matrix(data = NA, nrow = length(na_posicoes), ncol = ncol(agua)))

colnames(new_bacias) <- colnames(agua)

new_bacias$Group.1 <- rownames(estrutura)[na_posicoes]
new_bacias$X.urb <- urb[na_posicoes]

new_agua <- rbind(agua, new_bacias)

rownames(new_agua) <- new_agua$Group.1
urb_agua <- new_agua$X.urb

new_agua<- new_agua[,-c(1,2,3)]

#################################### GAM Imput #############################
library(mgcv)

predicted <- list()

families <- c("Gamma", "Gamma", "gaussian", "Gamma", "gaussian", "Gamma", "gaussian", "Gamma", "Gamma", "Gamma", "Gamma")
names(families) <- colnames(new_agua)

new_agua_imputed <- new_agua

for(i in 1:ncol(new_agua)){
  y <- new_agua[,i]
  gam_model <- gam(y ~ s(urb_agua), family = families[i])
  newdata <- data.frame(urb_agua = urb_agua[(length(urb_agua) - 3) :length(urb_agua) ])
  pred <- predict(gam_model, newdata = newdata, se = TRUE, type = "response")
  pred_fit <- pred$fit
  pred_se <- pred$se.fit
  set.seed(45)
  for(j in 1:length(pred_fit)){
    pred_fit[j] <- rnorm(1, mean = pred_fit[j], sd = pred_se[j])
  }
  predicted[[i]] <- pred_fit
  new_agua_imputed[(length(urb_agua) - 3) :length(urb_agua), i] <- pred_fit
}

##############################################################################

#Removendo bacias não usadas e reordenando a planilha de água de acordo com a planilha de estrutura

#pegando as posições na planilha água
posicoes <- match(rownames(estrutura), rownames(new_agua_imputed))

agua_selecionadas <- new_agua_imputed[posicoes,]

agua <- data.frame(bacia = rownames(estrutura), urb = urb, agua_selecionadas)

agua <- agua[,-c(1:2)]




#############################################################################
################################ DESCARGA ###################################
#############################################################################


#Adicionando descarga para planilha de estrutura e comparando com a outra descarga
estrutura$descarga_.L.s_sal <- agua$Discharge_.L.s.
agua <- agua[,colnames(agua) != "Discharge_.L.s."]

plot(x = estrutura$descarga_.L.s_sal, 
     y = estrutura$descarga_vel_dim,
     xlab = "Sal",
     ylab = "Calculado com velocidade")

#outliers
rownames(estrutura)[estrutura$descarga_.L.s_sal == max(na.omit(estrutura$descarga_.L.s_sal))]
#b637 #tem grande descarga de acordo com dados de sal, mas não tão grande de acordo com dados de dimensões e velocidade
rownames(estrutura)[estrutura$descarga_vel_dim == max(na.omit(estrutura$descarga_vel_dim))]
#b579 #tem grande descarga de acordo com dados de dimensões e velocidade, mas baixa de acordo com dados de sal

#Número de velocidades com Sal sem medidas
length(which(is.na(estrutura$descarga_.L.s_sal)))

#vizualinando
data.frame(bacia = rownames(estrutura),
           sal = estrutura$descarga_.L.s_sal,
           velocidade = estrutura$descarga_vel_dim)


#Estimando dados faltantes de descarga com Sal
#Substituindo Inf por zero
estrutura$descarga_vel_dim[estrutura$descarga_vel_dim == Inf] <- 0
estrutura$velocidade_media[estrutura$velocidade_media == Inf] <- 0

#removendo outliers
sal <- estrutura$descarga_.L.s_sal
vel <- estrutura$descarga_vel_dim
data_descarga <- data.frame(bacia = rownames(estrutura),
                            sal,
                            vel)

data_descarga <- data_descarga[data_descarga$bacia != "b637" &  data_descarga$bacia != "b579",]

regressao <- lm(sal ~ vel, data_descarga)
plot(y = data_descarga$sal, 
     x = data_descarga$vel,
     ylab = "Sal",
     xlab = "Calculado com velocidade")
abline(regressao)

newdata <- data.frame(vel = estrutura$descarga_vel_dim)
predicted <- predict(regressao, newdata = newdata)

valores_a_serem_adicionados <- which(is.na(estrutura$descarga_.L.s_sal))

nova_descarga_sal <- estrutura$descarga_.L.s_sal

nova_descarga_sal[valores_a_serem_adicionados] <- predicted[valores_a_serem_adicionados]

estrutura$nova_descarga_sal <- nova_descarga_sal


############################################################################

#Removendo bacias não usadas e reordenando a planilha de bacia de acordo com a planilha de estrutura

#pegando as posições na planilha bacia
posicoes <- match(rownames(estrutura), bacia$riacho)
bacia_selecionadas <- bacia[posicoes,]

posicoes2 <- match(rownames(estrutura), bacia2$bacia)
bacia_selecionadas2<- bacia2[posicoes2,]

bacia_selecionadas <- cbind(bacia_selecionadas, bacia_selecionadas2)

var_interesse_bacia <- c(#"bacia",
                        "Area_ha",
                        "FOR_2021",
                        "URB_2021",
                        "AGR_2021",
                        "Ic",
                        "Kc",
                        "Declividade_av",
                        "Altitude_av")

bacia <- bacia_selecionadas[,match(var_interesse_bacia, colnames(bacia_selecionadas))]
rownames(bacia) <- bacia_selecionadas2$bacia

nrow(bacia)

##############################################################################

#Verificando
nrow(agua)
nrow(bacia)
nrow(estrutura)
nrow(dados_peixes)

colnames(agua)
colnames(bacia)
colnames(estrutura)

data.frame(agua = rownames(agua),
           bacia = rownames(bacia),
           estrutura = rownames(estrutura),
           urb = names(urb))

##############################################################################

#Verificando de urbanização de delineamento é compatível com a informação da bacia. Acho que a do delineamento é de 2019 e da bacia de 2021
variaveis_urb_delinemaneto <- data.frame(urb_delinamento_2019 = urb,
                                         urb_bacia_2021 = bacia$URB_2021)

write.csv(variaveis_urb_delinemaneto, "urb_2019_2021.csv")

plot(x = urb, y = bacia$URB_2021, xlab = "Delineamento", ylab = "2021 Gabriel")
plot(x = urb, y = bacia_selecionadas$URB_2019, xlab = "Delineamento", ylab = "2019 Gabriel")
plot(x = urb, y = bacia_selecionadas$URB_2018, xlab = "Delineamento", ylab = "2018 Gabriel")
plot(x = urb, y = bacia_selecionadas$URB_2017, xlab = "Delineamento", ylab = "2017 Gabriel")


bacia_selecionadas

#Existe variação, adicionando urbanização delineamento pra planilha de bacia

bacia$urbano_delineamento <- urb

#############################################################################

#Juntando a planilha num grande planilhão

ncol_agua <- ncol(agua)
ncol_estrutura <- ncol(estrutura)
ncol_bacia <- ncol(bacia)

planilha_ambientais <- cbind(agua, bacia, estrutura)

write.csv(planilha_ambientais, "data/planilha_ambiental_assembleias.csv")
write.csv(agua, "data/planilha_agua_assembleias.csv")
write.csv(bacia, "data/planilha_bacia_assembleias.csv")
write.csv(estrutura, "data/planilha_estrutura_assembleias.csv")


#############################################################################


