estrutura <- read.csv("estrutura.csv")
urbanizacao <- read.csv("urbanização_floresta.csv")

## Separando planilhas com os dados de composição de econoto, substrato, etc.
comp_ecotono <- estrutura[-1,estrutura[1,] == "comp_ecotono"]
comp_ecotono[] <- lapply(comp_ecotono, as.numeric)

substrato <- estrutura[-1,estrutura[1,] == "substrato"]
substrato[] <- lapply(substrato, as.numeric)

tipo_de_canal <- estrutura[-1,estrutura[1,] == "tipo_de_canal"]
tipo_de_canal[] <- lapply(tipo_de_canal, as.numeric)

encaixe <- estrutura[-1,estrutura[1,] == "encaixe"]
encaixe[] <- lapply(encaixe, as.numeric)

pert_zona_riparia <- estrutura[-1,estrutura[1,] == "pert_zona_riparia"]
pert_zona_riparia[] <- lapply(pert_zona_riparia, as.numeric)

comp_zona_riparia <- estrutura[-1,estrutura[1,] == "comp_zona_riparia"]
comp_zona_riparia[] <- lapply(comp_zona_riparia, as.numeric)

estrutura_dentro_do_canal <- estrutura[-1,estrutura[1,] == "estrutura_dentro_do_canal"]
estrutura_dentro_do_canal[] <- lapply(estrutura_dentro_do_canal, as.numeric)


#Transformando dossel e sinuosidade em numerico
sinuosidade <- as.numeric(estrutura[-1,estrutura[1,] == "sinuosidade"])
dossel <- as.numeric(estrutura[-1,estrutura[1,] == "dossel"])
prof_media <- as.numeric(estrutura$prof_media[-1])
prof_cv <- as.numeric(estrutura$prof_cv[-1])
larg_media <- as.numeric(estrutura$larg_media[-1])
larg_cv <- as.numeric(estrutura$larg_cv[-1])
larg_media <- as.numeric(estrutura$larg_media[-1])
velocidade_media <- as.numeric(estrutura$velocidade_media[-1])
velocidade_cv <- as.numeric(estrutura$velocidade_cv[-1])
padrao_do_canal <- estrutura$padrao_canal[-1]
tipo_de_encaixe <- estrutura$tipo_de_encaixe[-1]







############################################################################
#Sumarizando trechos em variaveis unicas para cada riacho


bacia <-estrutura[-1,1]



#em alguns casos não faremos uma média, mas tiraremos a moda
getmode <- function(v) {
  v<- na.omit(v)
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

sinuosidade_summ <- tapply(sinuosidade, INDEX = bacia, getmode)

encaixe_summ <- aggregate(encaixe, by = list(bacia), FUN = getmode)

padrao_do_canal_summ <- aggregate(padrao_do_canal, by = list(bacia), FUN = getmode)

tipo_de_encaixe_summ <- aggregate(tipo_de_encaixe, by = list(bacia), FUN = getmode)



#médias
#univariados
dossel_summ <- aggregate(dossel, by = list(bacia), FUN = mean, na.rm = TRUE)
prof_media_summ <- aggregate(prof_media, by = list(bacia), FUN = mean, na.rm = TRUE)
prof_cv_summ <- aggregate(prof_cv, by = list(bacia), FUN = mean, na.rm = TRUE)
larg_media_summ <- aggregate(larg_media, by = list(bacia), FUN = mean, na.rm = TRUE)
larg_cv_summ <- aggregate(larg_cv, by = list(bacia), FUN = mean, na.rm = TRUE)
velocidade_media_summ <- aggregate(velocidade_media, by = list(bacia), FUN = mean, na.rm = TRUE)
velocidade_cv_summ <- aggregate(velocidade_cv, by = list(bacia), FUN = mean, na.rm = TRUE)

#multivariados
comp_ecotono_summ <- aggregate(comp_ecotono, by = list(bacia), FUN = mean, na.rm = TRUE)
rownames(comp_ecotono_summ) <- comp_ecotono_summ$Group.1

substrato_summ <- aggregate(substrato, by = list(bacia), FUN = mean, na.rm = TRUE)
rownames(substrato_summ) <- substrato_summ$Group.1

tipo_de_canal_summ <- aggregate(tipo_de_canal, by = list(bacia), FUN = mean, na.rm = TRUE)
rownames(tipo_de_canal_summ) <- tipo_de_canal_summ$Group.1

pert_zona_riparia_summ <- aggregate(pert_zona_riparia, by = list(bacia), FUN = mean, na.rm = TRUE)
rownames(pert_zona_riparia_summ) <- pert_zona_riparia_summ$Group.1

comp_zona_riparia_summ <- aggregate(comp_zona_riparia, by = list(bacia), FUN = mean, na.rm = TRUE)
rownames(comp_zona_riparia_summ) <- comp_zona_riparia_summ$Group.1

estrutura_dentro_do_canal_summ <- aggregate(estrutura_dentro_do_canal, by = list(bacia), FUN = mean, na.rm = TRUE)
rownames(estrutura_dentro_do_canal_summ) <- estrutura_dentro_do_canal_summ$Group.1



#Juntando tudo novamente
dados_estrutura_sumarizados <- data.frame(bacia = c(NA,estrutura_dentro_do_canal_summ$Group.1),
           urbanizacao = c(NA, urbanizacao$urbana),
           dossel = c(NA,dossel_summ[,-1]),
           prof_media = c(NA,prof_media_summ[,-1]),
           prof_cv = c(NA,prof_cv_summ[,-1]),
           larg_media = c(NA,larg_media_summ[,-1]), 
           larg_cv = c(NA,larg_cv_summ[,-1]),
           velocidade_media = c(NA,velocidade_media_summ[,-1]),
           velocidade_cv = c(NA, velocidade_cv_summ[,-1]),
           rbind(rep("comp_ecotono", ncol(comp_ecotono_summ[,-1])),comp_ecotono_summ[,-1]),
           rbind(rep("substrato", ncol(substrato_summ[,-1])),substrato_summ[,-1]),
           rbind(rep("tipo_de_canal", ncol(tipo_de_canal_summ[,-1])),tipo_de_canal_summ[,-1]),
           rbind(rep("pert_zona_riparia", ncol(pert_zona_riparia_summ[,-1])),pert_zona_riparia_summ[,-1]),
           rbind(rep("comp_zona_riparia", ncol(comp_zona_riparia_summ[,-1])),comp_zona_riparia_summ[,-1]),
           rbind(rep("estrutura_dentro_do_canal", ncol(estrutura_dentro_do_canal_summ[,-1])),estrutura_dentro_do_canal_summ[,-1]))


write.csv(dados_estrutura_sumarizados,"dados_estrutura_sumarizados.csv")
