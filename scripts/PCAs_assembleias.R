agua <-read.csv("data/planilha_agua_assembleias.csv", row.names = 1)
estrutura <-read.csv("data/planilha_estrutura_assembleias.csv", row.names = 1)
bacia <- read.csv("data/planilha_bacia_assembleias.csv", row.names = 1)


#removendo urbanização da bacia
urb <- bacia$urbano_delineamento
bacia <- bacia[colnames(bacia) != "urbano_delineamento" &
                 colnames(bacia) != "urbano_2021_Total"]
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




################################################################################################################
#PCA estrutura
pca_estrutura <- rda(estrutura_st)

importance_estrutura <- round(pca_estrutura$CA$eig/sum(pca_estrutura$CA$eig),2)
Eigenvalues_estrutura <- data.frame(autovalores = pca_estrutura$CA$eig,
                                    importance = importance_estrutura)

sum(importance_estrutura[1:5])
estrutura_PCs <- pca_estrutura$CA$u
estrutura_loadings <- pca_estrutura$CA$v

estrutura_loadings_filtrados <- estrutura_loadings[which(estrutura_loadings[,1] > 0.2 | estrutura_loadings[,1] < -0.2 |
                                                 estrutura_loadings[,2] > 0.2 | estrutura_loadings[,2] < -0.2),]


write.csv(Eigenvalues_estrutura, "data/pcas_amb/estrutura_autovalores.csv")
write.csv(estrutura_PCs, "data/pcas_amb/estrutura_PCs.csv")
write.csv(estrutura_loadings, "data/pcas_amb/estrutura_loadings.csv")




####################################################################################################
#Plot parameters
scaler <- min(max(abs(estrutura_PCs[, 1]))/max(abs(estrutura_loadings_filtrados[,1])),
              max(abs(estrutura_PCs[, 2]))/max(abs(estrutura_loadings_filtrados[,2])))

estrutura_loadings_filtrados_new <- estrutura_loadings_filtrados * scaler * 0.8

#estrutura_PCs <- jitter(estrutura_PCs, amount = 0.1)

pc1_label <- paste("PC1 (",round(importance_estrutura[1]*100,2),"%)",sep = "")
pc2_label <- paste("PC2 (",round(importance_estrutura[2]*100,2),"%)",sep = "")

xmin <- min(c(estrutura_PCs[,1], estrutura_loadings_filtrados_new[,1]))*1.1
xmax <- max(c(estrutura_PCs[,1], estrutura_loadings_filtrados_new[,1]))*1.1
ymin <- min(c(estrutura_PCs[,2], estrutura_loadings_filtrados_new[,2]))*1.1
ymax <- max(c(estrutura_PCs[,2], estrutura_loadings_filtrados_new[,2]))*1.1



#Plot########################################################################################
pdf("plots/pca_estrutura.pdf", height = 3.5, width = 3.5, pointsize = 5)
par(bg = "white", mar = c(4,4,2,0.1), bty = "o", cex = 1.25)
plot(estrutura_PCs[,1], estrutura_PCs[,2], xlim = c(xmin,xmax), ylim = c(ymin, ymax),
     type = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "")

abline(h = 0, v = 0, lty = 2)

library(scales)
library(shape)


pal <- col_numeric(palette = c("white", "black"), domain = urb, na.color = "grey50", alpha = FALSE, reverse = FALSE)
col <-pal(urb)

Arrows(x0 <- rep(0, nrow(estrutura_loadings_filtrados_new)),
       y0 <- rep(0, nrow(estrutura_loadings_filtrados_new)),
       x1 <- estrutura_loadings_filtrados_new[,1],
       y1 <- estrutura_loadings_filtrados_new[,2], arr.type = "triangle", arr.length = 0.4, col = "brown", lwd = 1.5)

points(estrutura_PCs[,1],estrutura_PCs[,2], col = "black", bg = col, pch = 21, cex = 1.5)
#text(estrutura_PCs[,1],estrutura_PCs[,2], labels = rownames(estrutura_PCs))

text(x = estrutura_loadings_filtrados_new[,1], y = estrutura_loadings_filtrados_new[,2], labels = rownames(estrutura_loadings_filtrados_new), cex = 0.8)

axis(1, cex.axis = 1.25)
axis(2, cex.axis = 1.25, las = 2)
title(xlab = pc1_label, cex.lab = 1.4, line = 2.75)
title(ylab = pc2_label, cex.lab = 1.4, line = 2.75)
title(main = "Estrutura", line = 0.5, adj = 0, cex.main = 1.75)
dev.off()

############################################################################################




###########################################################################################
###########################################################################################


################################################################################################################
#PCA bacia
library(vegan)
pca_bacia <- rda(bacia_st)

importance_bacia <- round(pca_bacia$CA$eig/sum(pca_bacia$CA$eig),2)
Eigenvalues_bacia <- data.frame(autovalores = pca_bacia$CA$eig,
                                    importance = importance_bacia)

sum(importance_bacia[1:3])
bacia_PCs <- pca_bacia$CA$u
bacia_loadings <- pca_bacia$CA$v

bacia_loadings_filtrados <- bacia_loadings[which(bacia_loadings[,1] > 0.2 | bacia_loadings[,1] < -0.2 |
                                                 bacia_loadings[,2] > 0.2 | bacia_loadings[,2] < -0.2),]


write.csv(Eigenvalues_bacia, "data/pcas_amb/bacia_autovalores.csv")
write.csv(bacia_PCs, "data/pcas_amb/bacia_PCs.csv")
write.csv(bacia_loadings, "data/pcas_amb/bacia_loadings.csv")




####################################################################################################
#Plot parameters
scaler <- min(max(abs(bacia_PCs[, 1]))/max(abs(bacia_loadings_filtrados[,1])),
              max(abs(bacia_PCs[, 2]))/max(abs(bacia_loadings_filtrados[,2])))

bacia_loadings_filtrados_new <- bacia_loadings_filtrados * scaler * 0.8

#bacia_PCs <- jitter(bacia_PCs, amount = 0.1)

pc1_label <- paste("PC1 (",round(importance_bacia[1]*100,2),"%)",sep = "")
pc2_label <- paste("PC2 (",round(importance_bacia[2]*100,2),"%)",sep = "")

xmin <- min(c(bacia_PCs[,1], bacia_loadings_filtrados_new[,1]))*1.1
xmax <- max(c(bacia_PCs[,1], bacia_loadings_filtrados_new[,1]))*1.1
ymin <- min(c(bacia_PCs[,2], bacia_loadings_filtrados_new[,2]))*1.1
ymax <- max(c(bacia_PCs[,2], bacia_loadings_filtrados_new[,2]))*1.1



#Plot########################################################################################
pdf("plots/pca_bacia.pdf", height = 3.5, width = 3.5, pointsize = 5)
par(bg = "white", mar = c(4,4,2,0.1), bty = "o", cex = 1.25)
plot(bacia_PCs[,1], bacia_PCs[,2], xlim = c(xmin,xmax), ylim = c(ymin, ymax),
     type = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "")

abline(h = 0, v = 0, lty = 2)

library(scales)
library(shape)


pal <- col_numeric(palette = c("white", "black"), domain = urb, na.color = "grey50", alpha = FALSE, reverse = FALSE)
col <-pal(urb)

Arrows(x0 <- rep(0, nrow(bacia_loadings_filtrados_new)),
       y0 <- rep(0, nrow(bacia_loadings_filtrados_new)),
       x1 <- bacia_loadings_filtrados_new[,1],
       y1 <- bacia_loadings_filtrados_new[,2], arr.type = "triangle", arr.length = 0.4, col = "brown", lwd = 1.5)

points(bacia_PCs[,1],bacia_PCs[,2], col = "black", bg = col, pch = 21, cex = 1.5)
#text(bacia_PCs[,1],bacia_PCs[,2], labels = rownames(bacia_PCs))

text(x = bacia_loadings_filtrados_new[,1], y = bacia_loadings_filtrados_new[,2], labels = rownames(bacia_loadings_filtrados_new), cex = 0.8)

axis(1, cex.axis = 1.25)
axis(2, cex.axis = 1.25, las = 2)
title(xlab = pc1_label, cex.lab = 1.4, line = 2.75)
title(ylab = pc2_label, cex.lab = 1.4, line = 2.75)
title(main = "Bacia", line = 0.5, adj = 0, cex.main = 1.75)
dev.off()

############################################################################################







###########################################################################################
###########################################################################################


################################################################################################################
#PCA agua
library(vegan)
pca_agua <- rda(agua_st)

importance_agua <- round(pca_agua$CA$eig/sum(pca_agua$CA$eig),2)
Eigenvalues_agua <- data.frame(autovalores = pca_agua$CA$eig,
                                importance = importance_agua)

sum(importance_agua[1:2])
agua_PCs <- pca_agua$CA$u
agua_loadings <- pca_agua$CA$v

agua_loadings_filtrados <- agua_loadings[which(agua_loadings[,1] > 0.2 | agua_loadings[,1] < -0.2 |
                                                   agua_loadings[,2] > 0.2 | agua_loadings[,2] < -0.2),]


write.csv(Eigenvalues_agua, "data/pcas_amb/agua_autovalores.csv")
write.csv(agua_PCs, "data/pcas_amb/agua_PCs.csv")
write.csv(agua_loadings, "data/pcas_amb/agua_loadings.csv")




####################################################################################################
#Plot parameters
scaler <- min(max(abs(agua_PCs[, 1]))/max(abs(agua_loadings_filtrados[,1])),
              max(abs(agua_PCs[, 2]))/max(abs(agua_loadings_filtrados[,2])))

agua_loadings_filtrados_new <- agua_loadings_filtrados * scaler * 0.8

#agua_PCs <- jitter(agua_PCs, amount = 0.1)

pc1_label <- paste("PC1 (",round(importance_agua[1]*100,2),"%)",sep = "")
pc2_label <- paste("PC2 (",round(importance_agua[2]*100,2),"%)",sep = "")

xmin <- min(c(agua_PCs[,1], agua_loadings_filtrados_new[,1]))*1.1
xmax <- max(c(agua_PCs[,1], agua_loadings_filtrados_new[,1]))*1.1
ymin <- min(c(agua_PCs[,2], agua_loadings_filtrados_new[,2]))*1.1
ymax <- max(c(agua_PCs[,2], agua_loadings_filtrados_new[,2]))*1.1



#Plot########################################################################################
pdf("plots/pca_agua.pdf", height = 3.5, width = 3.5, pointsize = 5)
par(bg = "white", mar = c(4,4,2,0.1), bty = "o", cex = 1.25)
plot(agua_PCs[,1], agua_PCs[,2], xlim = c(xmin,xmax), ylim = c(ymin, ymax),
     type = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "")

abline(h = 0, v = 0, lty = 2)

library(scales)
library(shape)


pal <- col_numeric(palette = c("white", "black"), domain = urb, na.color = "grey50", alpha = FALSE, reverse = FALSE)
col <-pal(urb)

Arrows(x0 <- rep(0, nrow(agua_loadings_filtrados_new)),
       y0 <- rep(0, nrow(agua_loadings_filtrados_new)),
       x1 <- agua_loadings_filtrados_new[,1],
       y1 <- agua_loadings_filtrados_new[,2], arr.type = "triangle", arr.length = 0.4, col = "brown", lwd = 1.5)

points(agua_PCs[,1],agua_PCs[,2], col = "black", bg = col, pch = 21, cex = 1.5)
#text(agua_PCs[,1],agua_PCs[,2], labels = rownames(agua_PCs))

text(x = agua_loadings_filtrados_new[,1], y = agua_loadings_filtrados_new[,2], labels = rownames(agua_loadings_filtrados_new), cex = 0.8)

axis(1, cex.axis = 1.25)
axis(2, cex.axis = 1.25, las = 2)
title(xlab = pc1_label, cex.lab = 1.4, line = 2.75)
title(ylab = pc2_label, cex.lab = 1.4, line = 2.75)
title(main = "Água", line = 0.5, adj = 0, cex.main = 1.75)
dev.off()

############################################################################################



