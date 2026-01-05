sinuosidade_summ
dossel_summ
prof_media_summ
prof_cv_summ
larg_media_summ
larg_cv_summ
velocidade_media_summ
velocidade_cv_summ
padrao_do_canal_summ
tipo_de_encaixe_summ

urbanizacao <- dados_estrutura_sumarizados$urbanizacao[-1]


#Sinuosidade
par(bg = "white", mar = c(4,4,2,0.1), bty = "o", cex = 1.25)
plot(sinuosidade_summ ~ urbanizacao, xlab = "", ylab = "", type = "n")
title(xlab = "Urbanização", cex.lab = 1.4, line = 2.75)
title(ylab = "Sinuosidade", cex.lab = 1.4, line = 2.75)
points(sinuosidade_summ ~ urbanizacao, col = "black", bg = "grey60", pch = 21, cex = 1.5)


#Dossel
par(bg = "white", mar = c(4,4,2,0.1), bty = "o", cex = 1.25)
plot(dossel_summ[,-1] ~ urbanizacao, xlab = "", ylab = "", type = "n")
title(xlab = "Urbanização", cex.lab = 1.4, line = 2.75)
title(ylab = "Dossel", cex.lab = 1.4, line = 2.75)
points(dossel_summ[,-1] ~ urbanizacao, col = "black", bg = "grey60", pch = 21, cex = 1.5)

#Profundidade
par(bg = "white", mar = c(4,4,2,0.1), bty = "o", cex = 1.25)
plot(prof_media_summ[,-1] ~ urbanizacao, xlab = "", ylab = "", type = "n")
title(xlab = "Urbanização", cex.lab = 1.4, line = 2.75)
title(ylab = "Profundiade Média (cm)", cex.lab = 1.4, line = 2.75)
points(prof_media_summ[,-1] ~ urbanizacao, col = "black", bg = "grey60", pch = 21, cex = 1.5)


#Largura
par(bg = "white", mar = c(4,4,2,0.1), bty = "o", cex = 1.25)
plot(larg_media_summ[,-1] ~ urbanizacao, xlab = "", ylab = "", type = "n")
title(xlab = "Urbanização", cex.lab = 1.4, line = 2.75)
title(ylab = "Largura Média (m)", cex.lab = 1.4, line = 2.75)
points(larg_media_summ[,-1] ~ urbanizacao, col = "black", bg = "grey60", pch = 21, cex = 1.5)


#Profundidade CV
par(bg = "white", mar = c(4,4,2,0.1), bty = "o", cex = 1.25)
plot(prof_cv_summ[,-1] ~ urbanizacao, xlab = "", ylab = "", type = "n")
title(xlab = "Urbanização", cex.lab = 1.4, line = 2.75)
title(ylab = "Variação de profundidade (%)", cex.lab = 1.4, line = 2.75)
points(prof_cv_summ[,-1] ~ urbanizacao, col = "black", bg = "grey60", pch = 21, cex = 1.5)

#Largura CV
par(bg = "white", mar = c(4,4,2,0.1), bty = "o", cex = 1.25)
plot(larg_cv_summ[,-1] ~ urbanizacao, xlab = "", ylab = "", type = "n")
title(xlab = "Urbanização", cex.lab = 1.4, line = 2.75)
title(ylab = "Variação de largura (%)", cex.lab = 1.4, line = 2.75)
points(larg_cv_summ[,-1] ~ urbanizacao, col = "black", bg = "grey60", pch = 21, cex = 1.5)


#Velocidade 
par(bg = "white", mar = c(4,4,2,0.1), bty = "o", cex = 1.25)
plot(velocidade_media_summ[,-1] ~ urbanizacao, xlab = "", ylab = "", type = "n")
title(xlab = "Urbanização", cex.lab = 1.4, line = 2.75)
title(ylab = "Velocidade média (s/m)", cex.lab = 1.4, line = 2.75)
points(velocidade_media_summ[,-1] ~ urbanizacao, col = "black", bg = "grey60", pch = 21, cex = 1.5)



#Velocidade CV
par(bg = "white", mar = c(4,4,2,0.1), bty = "o", cex = 1.25)
plot(velocidade_cv_summ[,-1] ~ urbanizacao, xlab = "", ylab = "", type = "n")
title(xlab = "Urbanização", cex.lab = 1.4, line = 2.75)
title(ylab = "Variação de velocidade(%)", cex.lab = 1.4, line = 2.75)
points(velocidade_cv_summ[,-1] ~ urbanizacao, col = "black", bg = "grey60", pch = 21, cex = 1.5)



#Velocidade CV
par(bg = "white", mar = c(4,4,2,0.1), bty = "o", cex = 1.25)
plot(velocidade_cv_summ[,-1] ~ urbanizacao, xlab = "", ylab = "", type = "n")
title(xlab = "Urbanização", cex.lab = 1.4, line = 2.75)
title(ylab = "Variação de velocidade(%)", cex.lab = 1.4, line = 2.75)
points(velocidade_cv_summ[,-1] ~ urbanizacao, col = "black", bg = "grey60", pch = 21, cex = 1.5)


tipo_de_encaixe_summ
