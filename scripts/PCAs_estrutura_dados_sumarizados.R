
library(vegan)


###########################################################################################
###########################################################################################

#PCA Composição Ecotono
comp_ecotono_summ_rda<- rda(comp_ecotono_summ[,-1])
importance <- round(comp_ecotono_summ_rda$CA$eig/sum(comp_ecotono_summ_rda$CA$eig),2)
comp_ecotono_PC <- comp_ecotono_summ_rda$CA$u[,1]
pcas <- comp_ecotono_summ_rda$CA$u
loadings <- comp_ecotono_summ_rda$CA$v

loadings_filtrados <- loadings[which(loadings[,1] > 0.2 | loadings[,1] < -0.2 | loadings[,2] > 0.2 | loadings[,2] < -0.2),]

#Plot parameters
scaler <- min(max(abs(pcas[, 1]))/max(abs(loadings[,1])),
              max(abs(pcas[, 2]))/max(abs(loadings[,2])))

loadings_new <- loadings_filtrados * scaler * 0.8

#pcas <- jitter(pcas, amount = 0.1)

pc1_label <- paste("PC1 (",round(importance[1]*100,2),"%)",sep = "")
pc2_label <- paste("PC2 (",round(importance[2]*100,2),"%)",sep = "")

xmin <- min(c(pcas[,1], loadings_new[,1]))*1.1
xmax <- max(c(pcas[,1], loadings_new[,1]))*1.1
ymin <- min(c(pcas[,2], loadings_new[,2]))*1.1
ymax <- max(c(pcas[,2], loadings_new[,2]))*1.1



#Plot#############################################################################
pdf("Plots/PCA_composicao_ecotono_pontos.pdf", height = 3.5, width = 3.5, pointsize = 5)
par(bg = "white", mar = c(4,4,2,0.1), bty = "o", cex = 1.25)
plot(pcas[,1], pcas[,2], xlim = c(xmin,xmax), ylim = c(ymin, ymax),
     type = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "")

abline(h = 0, v = 0, lty = 2)

library(scales)
library(shape)

dados_estrutura_sumarizados$urbanizacao

pal <- col_numeric(palette = c("white", "black"), domain = dados_estrutura_sumarizados$urbanizacao[-1], na.color = "grey50", alpha = FALSE, reverse = FALSE)
col <-pal(dados_estrutura_sumarizados$urbanizacao[-1])

Arrows(x0 <- rep(0, nrow(loadings_new)),
       y0 <- rep(0, nrow(loadings_new)),
       x1 <- loadings_new[,1],
       y1 <- loadings_new[,2], arr.type = "triangle", arr.length = 0.4, col = "brown", lwd = 1.5)

points(pcas[,1],pcas[,2], col = "black", bg = col, pch = 21, cex = 1.5)
#text(pcas[,1],pcas[,2], labels = rownames(pcas))

text(x = loadings_new[,1], y = loadings_new[,2], labels = rownames(loadings_new))

axis(1, cex.axis = 1.25)
axis(2, cex.axis = 1.25, las = 2)
title(xlab = pc1_label, cex.lab = 1.4, line = 2.75)
title(ylab = pc2_label, cex.lab = 1.4, line = 2.75)
title(main = "Composição do Ecotono", line = 0.5, adj = 0, cex.main = 1.75)
dev.off()


#Plot#############################################################################
pdf("Plots/PCA_composicao_ecotono_nomes.pdf", height = 3.5, width = 3.5, pointsize = 5)
par(bg = "white", mar = c(4,4,2,0.1), bty = "o", cex = 1.25)
plot(pcas[,1], pcas[,2], xlim = c(xmin,xmax), ylim = c(ymin, ymax),
     type = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "")

abline(h = 0, v = 0, lty = 2)

library(scales)
library(shape)

dados_estrutura_sumarizados$urbanizacao

pal <- col_numeric(palette = c("white", "black"), domain = dados_estrutura_sumarizados$urbanizacao[-1], na.color = "grey50", alpha = FALSE, reverse = FALSE)
col <-pal(dados_estrutura_sumarizados$urbanizacao[-1])


Arrows(x0 <- rep(0, nrow(loadings_new)),
       y0 <- rep(0, nrow(loadings_new)),
       x1 <- loadings_new[,1],
       y1 <- loadings_new[,2], arr.type = "triangle", arr.length = 0.4, col = "brown", lwd = 1.5)

#points(pcas[,1],pcas[,2], col = "black", bg = col, pch = 21, cex = 1.5)
text(pcas[,1],pcas[,2], labels = rownames(pcas))


text(x = loadings_new[,1], y = loadings_new[,2], labels = rownames(loadings_new))

axis(1, cex.axis = 1.25)
axis(2, cex.axis = 1.25, las = 2)
title(xlab = pc1_label, cex.lab = 1.4, line = 2.75)
title(ylab = pc2_label, cex.lab = 1.4, line = 2.75)
title(main = "Composição do Ecotono", line = 0.5, adj = 0, cex.main = 1.75)

dev.off()
#######################################################################################





###########################################################################################
###########################################################################################

#PCA Substrato
substrato_summ_rda<- rda(substrato_summ[,-1])
importance <- round(substrato_summ_rda$CA$eig/sum(substrato_summ_rda$CA$eig),2)
substrato_PC <- substrato_rda$CA$u[,1]
pcas <- substrato_rda$CA$u
loadings <- substrato_rda$CA$v

loadings_filtrados <- loadings[which(loadings[,1] > 0.2 | loadings[,1] < -0.2 | loadings[,2] > 0.2 | loadings[,2] < -0.2),]

#Plot parameters
scaler <- min(max(abs(pcas[, 1]))/max(abs(loadings[,1])),
              max(abs(pcas[, 2]))/max(abs(loadings[,2])))

loadings_new <- loadings_filtrados * scaler * 0.8

#pcas <- jitter(pcas, amount = 0.1)

pc1_label <- paste("PC1 (",round(importance[1]*100,2),"%)",sep = "")
pc2_label <- paste("PC2 (",round(importance[2]*100,2),"%)",sep = "")

xmin <- min(c(pcas[,1], loadings_new[,1]))*1.1
xmax <- max(c(pcas[,1], loadings_new[,1]))*1.1
ymin <- min(c(pcas[,2], loadings_new[,2]))*1.1
ymax <- max(c(pcas[,2], loadings_new[,2]))*1.1



#Plot#############################################################################
pdf("Plots/PCA_substrato_pontos.pdf", height = 3.5, width = 3.5, pointsize = 5)
par(bg = "white", mar = c(4,4,2,0.1), bty = "o", cex = 1.25)
plot(pcas[,1], pcas[,2], xlim = c(xmin,xmax), ylim = c(ymin, ymax),
     type = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "")

abline(h = 0, v = 0, lty = 2)

dados_estrutura_sumarizados$urbanizacao

pal <- col_numeric(palette = c("white", "black"), domain = dados_estrutura_sumarizados$urbanizacao[-1], na.color = "grey50", alpha = FALSE, reverse = FALSE)
col <-pal(dados_estrutura_sumarizados$urbanizacao[-1])

Arrows(x0 <- rep(0, nrow(loadings_new)),
       y0 <- rep(0, nrow(loadings_new)),
       x1 <- loadings_new[,1],
       y1 <- loadings_new[,2], arr.type = "triangle", arr.length = 0.4, col = "brown", lwd = 1.5)

points(pcas[,1],pcas[,2], col = "black", bg = col, pch = 21, cex = 1.5)
#text(pcas[,1],pcas[,2], labels = rownames(pcas))

text(x = loadings_new[,1], y = loadings_new[,2], labels = rownames(loadings_new))

axis(1, cex.axis = 1.25)
axis(2, cex.axis = 1.25, las = 2)
title(xlab = pc1_label, cex.lab = 1.4, line = 2.75)
title(ylab = pc2_label, cex.lab = 1.4, line = 2.75)
title(main = "Substrato", line = 0.5, adj = 0, cex.main = 1.75)
dev.off()


#Plot#############################################################################
pdf("Plots/PCA_substrato_nomes.pdf", height = 3.5, width = 3.5, pointsize = 5)
par(bg = "white", mar = c(4,4,2,0.1), bty = "o", cex = 1.25)
plot(pcas[,1], pcas[,2], xlim = c(xmin,xmax), ylim = c(ymin, ymax),
     type = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "")

abline(h = 0, v = 0, lty = 2)

library(scales)
library(shape)

dados_estrutura_sumarizados$urbanizacao

pal <- col_numeric(palette = c("white", "black"), domain = dados_estrutura_sumarizados$urbanizacao[-1], na.color = "grey50", alpha = FALSE, reverse = FALSE)
col <-pal(dados_estrutura_sumarizados$urbanizacao[-1])


Arrows(x0 <- rep(0, nrow(loadings_new)),
       y0 <- rep(0, nrow(loadings_new)),
       x1 <- loadings_new[,1],
       y1 <- loadings_new[,2], arr.type = "triangle", arr.length = 0.4, col = "brown", lwd = 1.5)

#points(pcas[,1],pcas[,2], col = "black", bg = col, pch = 21, cex = 1.5)
text(pcas[,1],pcas[,2], labels = rownames(pcas))


text(x = loadings_new[,1], y = loadings_new[,2], labels = rownames(loadings_new))

axis(1, cex.axis = 1.25)
axis(2, cex.axis = 1.25, las = 2)
title(xlab = pc1_label, cex.lab = 1.4, line = 2.75)
title(ylab = pc2_label, cex.lab = 1.4, line = 2.75)
title(main = "Substrato", line = 0.5, adj = 0, cex.main = 1.75)

dev.off()
####################################################################################










###########################################################################################
###########################################################################################

#PCA Substrato
substrato_rda<- rda(substrato_summ[,-1])
importance <- round(substrato_summ_rda$CA$eig/sum(substrato_summ_rda$CA$eig),2)
substrato_PC <- substrato_rda$CA$u[,1]
pcas <- substrato_rda$CA$u
loadings <- substrato_rda$CA$v

loadings_filtrados <- loadings[which(loadings[,1] > 0.2 | loadings[,1] < -0.2 | loadings[,2] > 0.2 | loadings[,2] < -0.2),]

#Plot parameters
scaler <- min(max(abs(pcas[, 1]))/max(abs(loadings[,1])),
              max(abs(pcas[, 2]))/max(abs(loadings[,2])))

loadings_new <- loadings_filtrados * scaler * 0.8

#pcas <- jitter(pcas, amount = 0.1)

pc1_label <- paste("PC1 (",round(importance[1]*100,2),"%)",sep = "")
pc2_label <- paste("PC2 (",round(importance[2]*100,2),"%)",sep = "")

xmin <- min(c(pcas[,1], loadings_new[,1]))*1.1
xmax <- max(c(pcas[,1], loadings_new[,1]))*1.1
ymin <- min(c(pcas[,2], loadings_new[,2]))*1.1
ymax <- max(c(pcas[,2], loadings_new[,2]))*1.1



#Plot#############################################################################
pdf("Plots/PCA_substrato_pontos.pdf", height = 3.5, width = 3.5, pointsize = 5)
par(bg = "white", mar = c(4,4,2,0.1), bty = "o", cex = 1.25)
plot(pcas[,1], pcas[,2], xlim = c(xmin,xmax), ylim = c(ymin, ymax),
     type = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "")

abline(h = 0, v = 0, lty = 2)

dados_estrutura_sumarizados$urbanizacao

pal <- col_numeric(palette = c("white", "black"), domain = dados_estrutura_sumarizados$urbanizacao[-1], na.color = "grey50", alpha = FALSE, reverse = FALSE)
col <-pal(dados_estrutura_sumarizados$urbanizacao[-1])

Arrows(x0 <- rep(0, nrow(loadings_new)),
       y0 <- rep(0, nrow(loadings_new)),
       x1 <- loadings_new[,1],
       y1 <- loadings_new[,2], arr.type = "triangle", arr.length = 0.4, col = "brown", lwd = 1.5)

points(pcas[,1],pcas[,2], col = "black", bg = col, pch = 21, cex = 1.5)
#text(pcas[,1],pcas[,2], labels = rownames(pcas))

text(x = loadings_new[,1], y = loadings_new[,2], labels = rownames(loadings_new))

axis(1, cex.axis = 1.25)
axis(2, cex.axis = 1.25, las = 2)
title(xlab = pc1_label, cex.lab = 1.4, line = 2.75)
title(ylab = pc2_label, cex.lab = 1.4, line = 2.75)
title(main = "Substrato", line = 0.5, adj = 0, cex.main = 1.75)
dev.off()


#Plot#############################################################################
pdf("Plots/PCA_substrato_nomes.pdf", height = 3.5, width = 3.5, pointsize = 5)
par(bg = "white", mar = c(4,4,2,0.1), bty = "o", cex = 1.25)
plot(pcas[,1], pcas[,2], xlim = c(xmin,xmax), ylim = c(ymin, ymax),
     type = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "")

abline(h = 0, v = 0, lty = 2)

library(scales)
library(shape)

dados_estrutura_sumarizados$urbanizacao

pal <- col_numeric(palette = c("white", "black"), domain = dados_estrutura_sumarizados$urbanizacao[-1], na.color = "grey50", alpha = FALSE, reverse = FALSE)
col <-pal(dados_estrutura_sumarizados$urbanizacao[-1])


Arrows(x0 <- rep(0, nrow(loadings_new)),
       y0 <- rep(0, nrow(loadings_new)),
       x1 <- loadings_new[,1],
       y1 <- loadings_new[,2], arr.type = "triangle", arr.length = 0.4, col = "brown", lwd = 1.5)

#points(pcas[,1],pcas[,2], col = "black", bg = col, pch = 21, cex = 1.5)
text(pcas[,1],pcas[,2], labels = rownames(pcas))


text(x = loadings_new[,1], y = loadings_new[,2], labels = rownames(loadings_new))

axis(1, cex.axis = 1.25)
axis(2, cex.axis = 1.25, las = 2)
title(xlab = pc1_label, cex.lab = 1.4, line = 2.75)
title(ylab = pc2_label, cex.lab = 1.4, line = 2.75)
title(main = "Substrato", line = 0.5, adj = 0, cex.main = 1.75)

dev.off()
####################################################################################






###########################################################################################
###########################################################################################

#PCA tipo_de_canal
tipo_de_canal_rda<- rda(tipo_de_canal_summ[,-1])
importance <- round(tipo_de_canal_rda$CA$eig/sum(tipo_de_canal_rda$CA$eig),2)
tipo_de_canal_PC <- tipo_de_canal_rda$CA$u[,1]
pcas <- tipo_de_canal_rda$CA$u
loadings <- tipo_de_canal_rda$CA$v

loadings_filtrados <- loadings[which(loadings[,1] > 0.2 | loadings[,1] < -0.2 | loadings[,2] > 0.2 | loadings[,2] < -0.2),]

#Plot parameters
scaler <- min(max(abs(pcas[, 1]))/max(abs(loadings[,1])),
              max(abs(pcas[, 2]))/max(abs(loadings[,2])))

loadings_new <- loadings_filtrados * scaler * 0.8

#pcas <- jitter(pcas, amount = 0.1)

pc1_label <- paste("PC1 (",round(importance[1]*100,2),"%)",sep = "")
pc2_label <- paste("PC2 (",round(importance[2]*100,2),"%)",sep = "")

xmin <- min(c(pcas[,1], loadings_new[,1]))*1.1
xmax <- max(c(pcas[,1], loadings_new[,1]))*1.1
ymin <- min(c(pcas[,2], loadings_new[,2]))*1.1
ymax <- max(c(pcas[,2], loadings_new[,2]))*1.1



#Plot#############################################################################
pdf("Plots/PCA_tipo_de_canal_pontos.pdf", height = 3.5, width = 3.5, pointsize = 5)
par(bg = "white", mar = c(4,4,2,0.1), bty = "o", cex = 1.25)
plot(pcas[,1], pcas[,2], xlim = c(xmin,xmax), ylim = c(ymin, ymax),
     type = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "")

abline(h = 0, v = 0, lty = 2)

dados_estrutura_sumarizados$urbanizacao

pal <- col_numeric(palette = c("white", "black"), domain = dados_estrutura_sumarizados$urbanizacao[-1], na.color = "grey50", alpha = FALSE, reverse = FALSE)
col <-pal(dados_estrutura_sumarizados$urbanizacao[-1])

Arrows(x0 <- rep(0, nrow(loadings_new)),
       y0 <- rep(0, nrow(loadings_new)),
       x1 <- loadings_new[,1],
       y1 <- loadings_new[,2], arr.type = "triangle", arr.length = 0.4, col = "brown", lwd = 1.5)

points(pcas[,1],pcas[,2], col = "black", bg = col, pch = 21, cex = 1.5)
#text(pcas[,1],pcas[,2], labels = rownames(pcas))

text(x = loadings_new[,1], y = loadings_new[,2], labels = rownames(loadings_new))

axis(1, cex.axis = 1.25)
axis(2, cex.axis = 1.25, las = 2)
title(xlab = pc1_label, cex.lab = 1.4, line = 2.75)
title(ylab = pc2_label, cex.lab = 1.4, line = 2.75)
title(main = "Tipo de canal", line = 0.5, adj = 0, cex.main = 1.75)
dev.off()


#Plot#############################################################################
pdf("Plots/PCA_tipo_de_canal_nomes.pdf", height = 3.5, width = 3.5, pointsize = 5)
par(bg = "white", mar = c(4,4,2,0.1), bty = "o", cex = 1.25)
plot(pcas[,1], pcas[,2], xlim = c(xmin,xmax), ylim = c(ymin, ymax),
     type = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "")

abline(h = 0, v = 0, lty = 2)

library(scales)
library(shape)

dados_estrutura_sumarizados$urbanizacao

pal <- col_numeric(palette = c("white", "black"), domain = dados_estrutura_sumarizados$urbanizacao[-1], na.color = "grey50", alpha = FALSE, reverse = FALSE)
col <-pal(dados_estrutura_sumarizados$urbanizacao[-1])


Arrows(x0 <- rep(0, nrow(loadings_new)),
       y0 <- rep(0, nrow(loadings_new)),
       x1 <- loadings_new[,1],
       y1 <- loadings_new[,2], arr.type = "triangle", arr.length = 0.4, col = "brown", lwd = 1.5)

#points(pcas[,1],pcas[,2], col = "black", bg = col, pch = 21, cex = 1.5)
text(pcas[,1],pcas[,2], labels = rownames(pcas))


text(x = loadings_new[,1], y = loadings_new[,2], labels = rownames(loadings_new))

axis(1, cex.axis = 1.25)
axis(2, cex.axis = 1.25, las = 2)
title(xlab = pc1_label, cex.lab = 1.4, line = 2.75)
title(ylab = pc2_label, cex.lab = 1.4, line = 2.75)
title(main = "Tiop de canal", line = 0.5, adj = 0, cex.main = 1.75)

dev.off()
####################################################################################





###########################################################################################
###########################################################################################

#PCA comp_zona_riparia
comp_zona_riparia_rda<- rda(comp_zona_riparia_summ[,-1])
importance <- round(comp_zona_riparia_rda$CA$eig/sum(comp_zona_riparia_rda$CA$eig),2)
comp_zona_riparia_PC <- comp_zona_riparia_rda$CA$u[,1]
pcas <- comp_zona_riparia_rda$CA$u
loadings <- comp_zona_riparia_rda$CA$v

loadings_filtrados <- loadings[which(loadings[,1] > 0.2 | loadings[,1] < -0.2 | loadings[,2] > 0.2 | loadings[,2] < -0.2),]

#Plot parameters
scaler <- min(max(abs(pcas[, 1]))/max(abs(loadings[,1])),
              max(abs(pcas[, 2]))/max(abs(loadings[,2])))

loadings_new <- loadings_filtrados * scaler * 0.8

#pcas <- jitter(pcas, amount = 0.1)

pc1_label <- paste("PC1 (",round(importance[1]*100,2),"%)",sep = "")
pc2_label <- paste("PC2 (",round(importance[2]*100,2),"%)",sep = "")

xmin <- min(c(pcas[,1], loadings_new[,1]))*1.1
xmax <- max(c(pcas[,1], loadings_new[,1]))*1.1
ymin <- min(c(pcas[,2], loadings_new[,2]))*1.1
ymax <- max(c(pcas[,2], loadings_new[,2]))*1.1



#Plot#############################################################################
pdf("Plots/PCA_comp_zona_riparia_pontos.pdf", height = 3.5, width = 3.5, pointsize = 5)
par(bg = "white", mar = c(4,4,2,0.1), bty = "o", cex = 1.25)
plot(pcas[,1], pcas[,2], xlim = c(xmin,xmax), ylim = c(ymin, ymax),
     type = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "")

abline(h = 0, v = 0, lty = 2)

dados_estrutura_sumarizados$urbanizacao

pal <- col_numeric(palette = c("white", "black"), domain = dados_estrutura_sumarizados$urbanizacao[-1], na.color = "grey50", alpha = FALSE, reverse = FALSE)
col <-pal(dados_estrutura_sumarizados$urbanizacao[-1])

Arrows(x0 <- rep(0, nrow(loadings_new)),
       y0 <- rep(0, nrow(loadings_new)),
       x1 <- loadings_new[,1],
       y1 <- loadings_new[,2], arr.type = "triangle", arr.length = 0.4, col = "brown", lwd = 1.5)

points(pcas[,1],pcas[,2], col = "black", bg = col, pch = 21, cex = 1.5)
#text(pcas[,1],pcas[,2], labels = rownames(pcas))

text(x = loadings_new[,1], y = loadings_new[,2], labels = rownames(loadings_new))

axis(1, cex.axis = 1.25)
axis(2, cex.axis = 1.25, las = 2)
title(xlab = pc1_label, cex.lab = 1.4, line = 2.75)
title(ylab = pc2_label, cex.lab = 1.4, line = 2.75)
title(main = "Composição da Zona Ripária", line = 0.5, adj = 0, cex.main = 1.75)
dev.off()


#Plot#############################################################################
pdf("Plots/PCA_comp_zona_riparia_nomes.pdf", height = 3.5, width = 3.5, pointsize = 5)
par(bg = "white", mar = c(4,4,2,0.1), bty = "o", cex = 1.25)
plot(pcas[,1], pcas[,2], xlim = c(xmin,xmax), ylim = c(ymin, ymax),
     type = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "")

abline(h = 0, v = 0, lty = 2)

library(scales)
library(shape)

dados_estrutura_sumarizados$urbanizacao

pal <- col_numeric(palette = c("white", "black"), domain = dados_estrutura_sumarizados$urbanizacao[-1], na.color = "grey50", alpha = FALSE, reverse = FALSE)
col <-pal(dados_estrutura_sumarizados$urbanizacao[-1])


Arrows(x0 <- rep(0, nrow(loadings_new)),
       y0 <- rep(0, nrow(loadings_new)),
       x1 <- loadings_new[,1],
       y1 <- loadings_new[,2], arr.type = "triangle", arr.length = 0.4, col = "brown", lwd = 1.5)

#points(pcas[,1],pcas[,2], col = "black", bg = col, pch = 21, cex = 1.5)
text(pcas[,1],pcas[,2], labels = rownames(pcas))


text(x = loadings_new[,1], y = loadings_new[,2], labels = rownames(loadings_new))

axis(1, cex.axis = 1.25)
axis(2, cex.axis = 1.25, las = 2)
title(xlab = pc1_label, cex.lab = 1.4, line = 2.75)
title(ylab = pc2_label, cex.lab = 1.4, line = 2.75)
title(main = "Composição da Zona Ripária", line = 0.5, adj = 0, cex.main = 1.75)

dev.off()
####################################################################################









###########################################################################################
###########################################################################################

#PCA pert_zona_riparia
pert_zona_riparia_rda<- rda(pert_zona_riparia_summ[,-1])
importance <- round(pert_zona_riparia_rda$CA$eig/sum(pert_zona_riparia_rda$CA$eig),2)
pert_zona_riparia_PC <- pert_zona_riparia_rda$CA$u[,1]
pcas <- pert_zona_riparia_rda$CA$u
loadings <- pert_zona_riparia_rda$CA$v

loadings_filtrados <- loadings[which(loadings[,1] > 0.2 | loadings[,1] < -0.2 | loadings[,2] > 0.2 | loadings[,2] < -0.2),]

#Plot parameters
scaler <- min(max(abs(pcas[, 1]))/max(abs(loadings[,1])),
              max(abs(pcas[, 2]))/max(abs(loadings[,2])))

loadings_new <- loadings_filtrados * scaler * 0.8

#pcas <- jitter(pcas, amount = 0.1)

pc1_label <- paste("PC1 (",round(importance[1]*100,2),"%)",sep = "")
pc2_label <- paste("PC2 (",round(importance[2]*100,2),"%)",sep = "")

xmin <- min(c(pcas[,1], loadings_new[,1]))*1.1
xmax <- max(c(pcas[,1], loadings_new[,1]))*1.1
ymin <- min(c(pcas[,2], loadings_new[,2]))*1.1
ymax <- max(c(pcas[,2], loadings_new[,2]))*1.1



#Plot#############################################################################
pdf("Plots/PCA_pert_zona_riparia_pontos.pdf", height = 3.5, width = 3.5, pointsize = 5)
par(bg = "white", mar = c(4,4,2,0.1), bty = "o", cex = 1.25)
plot(pcas[,1], pcas[,2], xlim = c(xmin,xmax), ylim = c(ymin, ymax),
     type = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "")

abline(h = 0, v = 0, lty = 2)

dados_estrutura_sumarizados$urbanizacao

pal <- col_numeric(palette = c("white", "black"), domain = dados_estrutura_sumarizados$urbanizacao[-1], na.color = "grey50", alpha = FALSE, reverse = FALSE)
col <-pal(dados_estrutura_sumarizados$urbanizacao[-1])

Arrows(x0 <- rep(0, nrow(loadings_new)),
       y0 <- rep(0, nrow(loadings_new)),
       x1 <- loadings_new[,1],
       y1 <- loadings_new[,2], arr.type = "triangle", arr.length = 0.4, col = "brown", lwd = 1.5)

points(pcas[,1],pcas[,2], col = "black", bg = col, pch = 21, cex = 1.5)
#text(pcas[,1],pcas[,2], labels = rownames(pcas))

text(x = loadings_new[,1], y = loadings_new[,2], labels = rownames(loadings_new))

axis(1, cex.axis = 1.25)
axis(2, cex.axis = 1.25, las = 2)
title(xlab = pc1_label, cex.lab = 1.4, line = 2.75)
title(ylab = pc2_label, cex.lab = 1.4, line = 2.75)
title(main = "Perturbações na Zona Ripária", line = 0.5, adj = 0, cex.main = 1.75)
dev.off()


#Plot#############################################################################
pdf("Plots/PCA_pert_zona_riparia_nomes.pdf", height = 3.5, width = 3.5, pointsize = 5)
par(bg = "white", mar = c(4,4,2,0.1), bty = "o", cex = 1.25)
plot(pcas[,1], pcas[,2], xlim = c(xmin,xmax), ylim = c(ymin, ymax),
     type = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "")

abline(h = 0, v = 0, lty = 2)

library(scales)
library(shape)

dados_estrutura_sumarizados$urbanizacao

pal <- col_numeric(palette = c("white", "black"), domain = dados_estrutura_sumarizados$urbanizacao[-1], na.color = "grey50", alpha = FALSE, reverse = FALSE)
col <-pal(dados_estrutura_sumarizados$urbanizacao[-1])


Arrows(x0 <- rep(0, nrow(loadings_new)),
       y0 <- rep(0, nrow(loadings_new)),
       x1 <- loadings_new[,1],
       y1 <- loadings_new[,2], arr.type = "triangle", arr.length = 0.4, col = "brown", lwd = 1.5)

#points(pcas[,1],pcas[,2], col = "black", bg = col, pch = 21, cex = 1.5)
text(pcas[,1],pcas[,2], labels = rownames(pcas))


text(x = loadings_new[,1], y = loadings_new[,2], labels = rownames(loadings_new))

axis(1, cex.axis = 1.25)
axis(2, cex.axis = 1.25, las = 2)
title(xlab = pc1_label, cex.lab = 1.4, line = 2.75)
title(ylab = pc2_label, cex.lab = 1.4, line = 2.75)
title(main = "Perturbações na Zona Ripária", line = 0.5, adj = 0, cex.main = 1.75)

dev.off()
####################################################################################








###########################################################################################
###########################################################################################

#PCA estrutura_dentro_do_canal
estrutura_dentro_do_canal_rda<- rda(estrutura_dentro_do_canal_summ[,-1])
importance <- round(estrutura_dentro_do_canal_rda$CA$eig/sum(estrutura_dentro_do_canal_rda$CA$eig),2)
estrutura_dentro_do_canal_PC <- estrutura_dentro_do_canal_rda$CA$u[,1]
pcas <- estrutura_dentro_do_canal_rda$CA$u
loadings <- estrutura_dentro_do_canal_rda$CA$v

loadings_filtrados <- loadings[which(loadings[,1] > 0.2 | loadings[,1] < -0.2 | loadings[,2] > 0.2 | loadings[,2] < -0.2),]

#Plot parameters
scaler <- min(max(abs(pcas[, 1]))/max(abs(loadings[,1])),
              max(abs(pcas[, 2]))/max(abs(loadings[,2])))

loadings_new <- loadings_filtrados * scaler * 0.8

#pcas <- jitter(pcas, amount = 0.1)

pc1_label <- paste("PC1 (",round(importance[1]*100,2),"%)",sep = "")
pc2_label <- paste("PC2 (",round(importance[2]*100,2),"%)",sep = "")

xmin <- min(c(pcas[,1], loadings_new[,1]))*1.1
xmax <- max(c(pcas[,1], loadings_new[,1]))*1.1
ymin <- min(c(pcas[,2], loadings_new[,2]))*1.1
ymax <- max(c(pcas[,2], loadings_new[,2]))*1.1




#Plot#############################################################################
pdf("Plots/PCA_estrutura_dentro_do_canal_pontos.pdf", height = 3.5, width = 3.5, pointsize = 5)
par(bg = "white", mar = c(4,4,2,0.1), bty = "o", cex = 1.25)
plot(pcas[,1], pcas[,2], xlim = c(xmin,xmax), ylim = c(ymin, ymax),
     type = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "")

abline(h = 0, v = 0, lty = 2)

dados_estrutura_sumarizados$urbanizacao

pal <- col_numeric(palette = c("white", "black"), domain = dados_estrutura_sumarizados$urbanizacao[-1], na.color = "grey50", alpha = FALSE, reverse = FALSE)
col <-pal(dados_estrutura_sumarizados$urbanizacao[-1])

Arrows(x0 <- rep(0, nrow(loadings_new)),
       y0 <- rep(0, nrow(loadings_new)),
       x1 <- loadings_new[,1],
       y1 <- loadings_new[,2], arr.type = "triangle", arr.length = 0.4, col = "brown", lwd = 1.5)

points(pcas[,1],pcas[,2], col = "black", bg = col, pch = 21, cex = 1.5)
#text(pcas[,1],pcas[,2], labels = rownames(pcas))

text(x = loadings_new[,1], y = loadings_new[,2], labels = rownames(loadings_new))

axis(1, cex.axis = 1.25)
axis(2, cex.axis = 1.25, las = 2)
title(xlab = pc1_label, cex.lab = 1.4, line = 2.75)
title(ylab = pc2_label, cex.lab = 1.4, line = 2.75)
title(main = "Estruturas dentro do canal", line = 0.5, adj = 0, cex.main = 1.75)
dev.off()


#Plot#############################################################################
pdf("Plots/PCA_estrutura_dentro_do_canal_nomes.pdf", height = 3.5, width = 3.5, pointsize = 5)
par(bg = "white", mar = c(4,4,2,0.1), bty = "o", cex = 1.25)
plot(pcas[,1], pcas[,2], xlim = c(xmin,xmax), ylim = c(ymin, ymax),
     type = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "")

abline(h = 0, v = 0, lty = 2)

library(scales)
library(shape)

dados_estrutura_sumarizados$urbanizacao

pal <- col_numeric(palette = c("white", "black"), domain = dados_estrutura_sumarizados$urbanizacao[-1], na.color = "grey50", alpha = FALSE, reverse = FALSE)
col <-pal(dados_estrutura_sumarizados$urbanizacao[-1])


Arrows(x0 <- rep(0, nrow(loadings_new)),
       y0 <- rep(0, nrow(loadings_new)),
       x1 <- loadings_new[,1],
       y1 <- loadings_new[,2], arr.type = "triangle", arr.length = 0.4, col = "brown", lwd = 1.5)

#points(pcas[,1],pcas[,2], col = "black", bg = col, pch = 21, cex = 1.5)
text(pcas[,1],pcas[,2], labels = rownames(pcas))


text(x = loadings_new[,1], y = loadings_new[,2], labels = rownames(loadings_new))

axis(1, cex.axis = 1.25)
axis(2, cex.axis = 1.25, las = 2)
title(xlab = pc1_label, cex.lab = 1.4, line = 2.75)
title(ylab = pc2_label, cex.lab = 1.4, line = 2.75)
title(main = "Estruturas dentro do canal", line = 0.5, adj = 0, cex.main = 1.75)

dev.off()
####################################################################################






