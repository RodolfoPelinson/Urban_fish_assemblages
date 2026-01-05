functional_raw <- read.csv("data/functional_data_raw.csv")
functional_indices <- read.csv("data/functional_data.csv", row.names = 1)

Poecilia_vivipara <- rep(NA, ncol(functional_indices))
names(Poecilia_vivipara) <- colnames(functional_indices)
Poecilia_vivipara <- as.data.frame(t(data.frame(Poecilia_vivipara)))

str(functional_raw)
functional_raw[,3:ncol(functional_raw)] <- apply(functional_raw[3:ncol(functional_raw)], 2, as.numeric)
str(functional_raw)


Poecilia_vivipara_AVG <- colMeans(functional_raw[functional_raw$Especie == "Poecilia_vivipara",3:ncol(functional_raw)])

Poecilia_vivipara_AVG <- as.data.frame(rbind(Poecilia_vivipara_AVG))

Poecilia_vivipara$Compression.index <- Poecilia_vivipara_AVG$Altura_do_corpo_MBH/Poecilia_vivipara_AVG$Largura_do_corpo_MBW

Poecilia_vivipara$Depression.index <- Poecilia_vivipara_AVG$Altura_da_linha_media_BMH / Poecilia_vivipara_AVG$Altura_do_corpo_MBH

Poecilia_vivipara$Relative.depth <- Poecilia_vivipara_AVG$Altura_do_corpo_MBH / Poecilia_vivipara_AVG$Comprimento_padrao_SL

Poecilia_vivipara$Fineness.ratio <- Poecilia_vivipara_AVG$Comprimento_padrao_SL / sqrt(Poecilia_vivipara_AVG$Altura_do_corpo_MBH * Poecilia_vivipara_AVG$Largura_do_corpo_MBW)

Poecilia_vivipara$Relative.length.of.caudal.peduncle <- Poecilia_vivipara_AVG$Comp_pedun_caudal_CPdL / Poecilia_vivipara_AVG$Comprimento_padrao_SL

Poecilia_vivipara$Relative.height.of.caudal.peduncle <- Poecilia_vivipara_AVG$Altura_do_pedun_caudal_CPdH / Poecilia_vivipara_AVG$Altura_do_corpo_MBH

Poecilia_vivipara$Relative.width.of.caudal.peduncle <- Poecilia_vivipara_AVG$Largura_do_pedun_caudal_CPdW / Poecilia_vivipara_AVG$Largura_do_corpo_MBW

Poecilia_vivipara$Relative.area.of.dorsal.fin <- Poecilia_vivipara_AVG$Area_da_nadadeira_dorsal_DA / (Poecilia_vivipara_AVG$Comprimento_padrao_SL^2)

Poecilia_vivipara$Relative.area.of.caudal.fin <- Poecilia_vivipara_AVG$Area_da_nadadeira_caudal_CA / (Poecilia_vivipara_AVG$Comprimento_padrao_SL^2)

Poecilia_vivipara$Relative.area.of.pectoral.fin <- Poecilia_vivipara_AVG$Area_da_nadadeira_peitoral_PtA / (Poecilia_vivipara_AVG$Comprimento_padrao_SL^2)

Poecilia_vivipara$Aspect.ratio.of.pectoral.fin <- (Poecilia_vivipara_AVG$Comp_da_peitoral_PtL^2) / Poecilia_vivipara_AVG$Area_da_nadadeira_peitoral_PtA

Poecilia_vivipara$Relative.length.of.head <- Poecilia_vivipara_AVG$Comp_da_cabeca_HdL / Poecilia_vivipara_AVG$Comprimento_padrao_SL

Poecilia_vivipara$Relative.height.of.head <- Poecilia_vivipara_AVG$Altura_da_cabeca_HdH / Poecilia_vivipara_AVG$Altura_do_corpo_MBH

Poecilia_vivipara$Relative.width.of.mouth <- Poecilia_vivipara_AVG$Largura_da_boca_MW / Poecilia_vivipara_AVG$Largura_do_corpo_MBW

Poecilia_vivipara$Eye.position <- Poecilia_vivipara_AVG$Altura_da_linha_media_do_olho_EH / Poecilia_vivipara_AVG$Altura_da_cabeca_HdH

Poecilia_vivipara$Relative.height.of.mouth <- Poecilia_vivipara_AVG$Altura_da_boca_MH / Poecilia_vivipara_AVG$Altura_do_corpo_MBH

Poecilia_vivipara$Relative.width.of.head <- Poecilia_vivipara_AVG$Largura_da_cabeca_HdW / Poecilia_vivipara_AVG$Largura_do_corpo_MBW

Poecilia_vivipara$Relative.area.of.eye <- Poecilia_vivipara_AVG$Area_do_olho_EA / (Poecilia_vivipara_AVG$Comprimento_padrao_SL^2)

Poecilia_vivipara$Aspect.ratio.of.caudal.fin <- (Poecilia_vivipara_AVG$Altura_da_caudal_CH^2) / Poecilia_vivipara_AVG$Area_da_nadadeira_caudal_CA

Poecilia_vivipara$Relative.area.of.anal.fin <- Poecilia_vivipara_AVG$Area_da_nadadeira_anal_AA / (Poecilia_vivipara_AVG$Comprimento_padrao_SL^2)

Poecilia_vivipara$Aspect.ratio.of.anal.fin <- (Poecilia_vivipara_AVG$Comp_da_anal_AL^2) / Poecilia_vivipara_AVG$Area_da_nadadeira_anal_AA

Poecilia_vivipara$Relative.area.of.pectoral.fin <- Poecilia_vivipara_AVG$Area_da_nadadeira_peitoral_PtA / (Poecilia_vivipara_AVG$Comprimento_padrao_SL^2)

Poecilia_vivipara$Relative.area.of.pelvic.fin <- Poecilia_vivipara_AVG$Area_da_nadadeira_pelvica_PvA / (Poecilia_vivipara_AVG$Comprimento_padrao_SL^2)


Poecilia_vivipara

functional_indices <- rbind(functional_indices, Poecilia_vivipara)

write.csv(functional_indices, "data/functional_data.csv")

