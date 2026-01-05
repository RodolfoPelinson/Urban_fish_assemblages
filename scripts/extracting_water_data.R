water_quality_seca <- read.csv("data/water_quality_data.csv")

###################################################################
var_interesse_agua <- c("X.urb",
                        "chlorophyll_a",
                        "phycocyanin",
                        "Temperature_.oC.",
                        "DO_.mg.L.",
                        #"Conductivity_.uS.cm.",
                        "pH",
                        "turbidity_.NTU.",
                        "redox_potential_.mV.",
                        #"TOC",
                        "TC",
                        #"IC",
                        "TN",
                        #"TP_.Valderrama.",
                        "Discharge_.L.s.",
                        "SPC_.uS.cm.")

#Ajeitando planilhas variaveis da água

water_quality_agua <- water_quality_seca[,match(var_interesse_agua, colnames(water_quality_seca))]
water_quality_agua <- apply(water_quality_agua, 2, as.numeric)

water_quality_agua_mean <- aggregate(water_quality_agua, by = list(water_quality_seca$Catchment),  FUN = mean, na.rm = TRUE)
water_quality_agua_sd <- aggregate(water_quality_agua, by = list(water_quality_seca$Catchment),  FUN = sd, na.rm = TRUE)

for(i in 3:ncol(water_quality_agua_sd)){
  water_quality_agua_sd[,i][water_quality_agua_sd[,i] == 0] <- Inf
  water_quality_agua_sd[,i][water_quality_agua_sd[,i] == Inf] <- min(water_quality_agua_sd[,i], na.rm = TRUE)
  water_quality_agua_sd[,i][is.na(water_quality_agua_sd[,i])] <- mean(water_quality_agua_sd[,i], na.rm = TRUE)
  weights_agua <- water_quality_agua_sd
  weights_agua[,i] <- 1/weights_agua[,i]
}

write.csv(water_quality_agua_mean, "data/dados_ambientais_agua_2021.csv")
