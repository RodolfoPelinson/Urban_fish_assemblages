agua <-read.csv("data/planilha_agua_assembleias.csv", row.names = 1)
estrutura <-read.csv("data/planilha_estrutura_assembleias.csv", row.names = 1)
bacia <- read.csv("data/planilha_bacia_assembleias.csv", row.names = 1)
delineamento<- read.csv("data/delineamento.csv")


estrutura2 <- read.csv("data/dados_ambientais_estrutura_por_riacho.csv")
estrutura2 <- estrutura2[-1,]

dim(estrutura2)
dim(estrutura)


#removendo urbanização da bacia
urb <- bacia$urbano_delineamento
bacia <- bacia[colnames(bacia) != "urbano_delineamento" &
                 colnames(bacia) != "URB_2021"]
names(urb) <- rownames(bacia)


#removendo descargas extras da estrutura
estrutura <- estrutura[colnames(estrutura) != "descarga_vel_dim" &
                         colnames(estrutura) != "descarga_.L.s_sal"&
                         colnames(estrutura) != "descarga_.L.s_sal"]



#Padronizando variaveis
library(vegan)

estrutura <- estrutura[,colnames(estrutura) != "comp_zona_riparia_plantacao"] #Variavel só com zeros

agua_st <- decostand(agua, method = "stand")
estrutura_st <- decostand(estrutura, method = "stand")
bacia_st <- decostand(bacia, method = "stand")


#numero de colunas de cada uma
ncol_agua <- ncol(agua_st)
ncol_estrutura <- ncol(estrutura_st)
ncol_bacia <- ncol(bacia_st)



#############################################################
data.frame(rownames(estrutura), delineamento$bacia_id, estrutura2$bacia)

urb_group <- as.factor(delineamento$Grupo)


#canopy cover
dossel_mean <- tapply(estrutura$dossel, urb_group, mean)
dossel_sd <-tapply(estrutura$dossel, urb_group, sd)



#riparian vegetation (Arborea + arbustiva)
veg_arbust_arboreal_riparian <- estrutura$comp_zona_riparia_veg_arborea + estrutura$comp_zona_riparia_veg_arborea
veg_arbust_arboreal_riparian <- (veg_arbust_arboreal_riparian/max(veg_arbust_arboreal_riparian))*100
veg_arbust_arbor_mean <- tapply(veg_arbust_arboreal_riparian, urb_group, mean)
veg_arbust_arbor_sd <- tapply(veg_arbust_arboreal_riparian, urb_group, sd)

#channel simplification
canal_encaixado_canalizado_mean <- tapply(estrutura$tipo_de_encaixe_encaixado_canalizado*100, urb_group, mean)
canal_encaixado_canalizado_sd <- tapply(estrutura$tipo_de_encaixe_encaixado_canalizado*100, urb_group, sd)

#reduced sinuosity
sinuosidade_mean <- tapply(estrutura$sinuosidade, urb_group, mean)
sinuosidade_sd <- tapply(estrutura$sinuosidade, urb_group, sd)

#diminished mesohabitat diversity****
mesohabitat <- data.frame(corredeira = estrutura$tipo_de_canal_corredeira,
                          cascata = estrutura$tipo_de_canal_cascata,
                          fluxo_continuo = estrutura$tipo_de_canal_fluxo_continuo,
                          poco = estrutura$tipo_de_canal_poco)

library(vegan)
diversity_mesohabitat<- diversity(mesohabitat, index = "shannon")

diversity_mesohabitat_mean <- tapply(diversity_mesohabitat, urb_group, mean)
diversity_mesohabitat_sd <- tapply(diversity_mesohabitat, urb_group, sd)

fluxo_continuo_mean <- tapply(estrutura$tipo_de_canal_fluxo_continuo, urb_group, mean)
fluxo_continuo_sd <- tapply(estrutura$tipo_de_canal_fluxo_continuo, urb_group, sd)

corredeira_mean <- tapply(estrutura$tipo_de_canal_corredeira, urb_group, mean)
corredeira_sd <- tapply(estrutura$tipo_de_canal_corredeira, urb_group, sd)

poco_mean <- tapply(estrutura$tipo_de_canal_poco, urb_group, mean)
poco_sd <- tapply(estrutura$tipo_de_canal_poco, urb_group, sd)

cascata_mean <- tapply(estrutura$tipo_de_canal_cascata, urb_group, mean)
cascata_sd <- tapply(estrutura$tipo_de_canal_cascata, urb_group, sd)



#lower accumulation of woody debris and leaf litter
madeiras_e_bancos_de_folhas <- c(estrutura$estrutura_dentro_do_canal_pedaco_de_madeira_grande +
                                   estrutura$estrutura_dentro_do_canal_banco_de_folhas + 
                                   estrutura$estrutura_dentro_do_canal_pedaco_de_madeira_pequeno)
  

madeiras_e_bancos_de_folhas_mean <- tapply(madeiras_e_bancos_de_folhas, urb_group, mean)
madeiras_e_bancos_de_folhas_sd <- tapply(madeiras_e_bancos_de_folhas, urb_group, sd)

#temperature
temperature_mean <- tapply(agua$Temperature_.oC., urb_group, mean)
temperature_sd <- tapply(agua$Temperature_.oC., urb_group, sd)


#pH
pH_mean <- tapply(agua$pH, urb_group, mean)
pH_sd <- tapply(agua$pH, urb_group, sd)

#conductivity
SPC_mean <- tapply(agua$SPC_.uS.cm., urb_group, mean)
SPC_sd <- tapply(agua$SPC_.uS.cm., urb_group, sd)

#turbidity
turbidity_mean <- tapply(agua$turbidity_.NTU., urb_group, mean)
turbidity_sd <- tapply(agua$turbidity_.NTU., urb_group, sd)

#total dissolved carbon
TC_mean <- tapply(agua$TC, urb_group, mean)
TC_sd <- tapply(agua$TC, urb_group, sd)

#total dissolved nitrogen
TN_mean <- tapply(agua$TN, urb_group, mean)
TN_sd <- tapply(agua$TN, urb_group, sd)

#chlorophyll a
chlorophyll_a_mean <- tapply(agua$chlorophyll_a, urb_group, mean)
chlorophyll_a_sd <- tapply(agua$chlorophyll_a, urb_group, sd)

#dissolved oxygen 
DO_mean <- tapply(agua$DO_.mg.L., urb_group, mean)
DO_sd <- tapply(agua$DO_.mg.L., urb_group, sd)






means_per_group <- data.frame(dossel_mean,
           dossel_sd,
           
           veg_arbust_arbor_mean,
           veg_arbust_arbor_sd,
           
           canal_encaixado_canalizado_mean,
           canal_encaixado_canalizado_sd,
           
           sinuosidade_mean,
           sinuosidade_sd,
           
           diversity_mesohabitat_mean,
           diversity_mesohabitat_sd,
           
           fluxo_continuo_mean,
           fluxo_continuo_sd,
           
           corredeira_mean,
           corredeira_sd,
           
           poco_mean,
           poco_sd,
           
           cascata_mean,
           cascata_sd,
           
           madeiras_e_bancos_de_folhas_mean,
           madeiras_e_bancos_de_folhas_sd,
           
           temperature_mean,
           temperature_sd,
           
           pH_mean,
           pH_sd,
           
           SPC_mean,
           SPC_sd,
           
           turbidity_mean,
           turbidity_sd,
           
           TC_mean,
           TC_sd,
           
           TN_mean,
           TN_sd,
           
           chlorophyll_a_mean,
           chlorophyll_a_sd,
           
           DO_mean,
           DO_sd)

write.csv(means_per_group, "means_per_group_LUIS_nature_city.csv")
