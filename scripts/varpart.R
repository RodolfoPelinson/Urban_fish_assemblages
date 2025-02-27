
###############################
Eigenvalues_estrutura
sum(importance_estrutura[1:10])
estrutura_PCs <- pca_estrutura$CA$u
estrutura_loadings <- pca_estrutura$CA$v

Eigenvalues_bacia
sum(importance_bacia[1:3])
bacia_PCs <- pca_bacia$CA$u
bacia_loadings <- pca_bacia$CA$v

Eigenvalues_agua
sum(importance_agua[1])
agua_PCs <- pca_agua$CA$u
agua_loadings <- pca_agua$CA$v

###############################

assembleia_peixes <- dados_peixes[match(rownames(estrutura_PCs), rownames(dados_peixes),),]

nrow(estrutura_PCs)
nrow(agua_PCs)
nrow(bacia_PCs)
nrow(assembleia_peixes)

dist_assembleias <- vegdist(assembleia_peixes, method = "euclidean")

varp1 <- varpart(assembleia_peixes, estrutura_PCs[,1:2], agua_PCs[,1:2], bacia_PCs[,1:2], urb)
varp2 <- varpart(assembleia_peixes, estrutura_PCs[,1:2], agua_PCs[,1:2], bacia_PCs[,1:2], urb, transfo = "hel")
varp3 <- varpart(assembleia_peixes, estrutura_PCs[,1:10], agua_PCs[,1], bacia_PCs[,1:3], urb, transfo = "hel")




#########################
estrutura_rda <- rda(X = assembleia_peixes, Y =  estrutura_PCs[,1:2])
anova.cca(estrutura_rda, permutations = 10000)
RsquareAdj(estrutura_rda)

agua_rda <- rda(X = assembleia_peixes, Y =  agua_PCs[,1:2])
anova.cca(agua_rda, permutations = 10000)
RsquareAdj(agua_rda)

bacia_rda <- rda(X = assembleia_peixes, Y =  bacia_PCs[,1:2])
anova.cca(bacia_rda, permutations = 10000)
RsquareAdj(bacia_rda)

urb_rda <- rda(X = assembleia_peixes, Y =  urb)
anova.cca(urb_rda, permutations = 10000)
RsquareAdj(urb_rda)

#########################
only_estrutura_rda <- rda(X = assembleia_peixes, Y =  estrutura_PCs[,1:2], Z = data.frame(agua_PCs[,1:2], bacia_PCs[,1:2], urb))
anova.cca(only_estrutura_rda, permutations = 10000)
RsquareAdj(only_estrutura_rda)

only_agua_rda <- rda(X = assembleia_peixes, Y =  agua_PCs[,1:2], Z = data.frame(estrutura_PCs[,1:2], bacia_PCs[,1:2], urb))
anova.cca(only_agua_rda, permutations = 10000)
RsquareAdj(only_agua_rda)

only_bacia_rda <- rda(X = assembleia_peixes, Y =  bacia_PCs[,1:2], Z = data.frame(agua_PCs[,1:2], estrutura_PCs[,1:2], urb))
anova.cca(only_bacia_rda, permutations = 10000)
RsquareAdj(only_bacia_rda)

only_urb_rda <- rda(X = assembleia_peixes, Y =  urb, Z = data.frame(agua_PCs[,1:2], bacia_PCs[,1:2], estrutura_PCs[,1:2]))
anova.cca(only_urb_rda, permutations = 10000)
RsquareAdj(only_urb_rda)



pdf("plots/varpart.pdf", height = 3.5, width = 5, pointsize = 5)
plot(varp1, Xnames = c(paste("Estrutura", round(RsquareAdj(estrutura_rda)$adj.r.squared,2)),
                       paste("Água",round(RsquareAdj(agua_rda)$adj.r.squared,2)),
                       paste("Bacia",round(RsquareAdj(bacia_rda)$adj.r.squared,2)),
                       paste("Urbanização",round(RsquareAdj(urb_rda)$adj.r.squared,2))))
dev.off()











assembleia_peixes_hell <- decostand(assembleia_peixes, method = "hell")



#########################
estrutura_rda <- rda(X = assembleia_peixes_hell, Y =  estrutura_PCs[,1:2])
anova.cca(estrutura_rda, permutations = 10000)
RsquareAdj(estrutura_rda)

agua_rda <- rda(X = assembleia_peixes_hell, Y =  agua_PCs[,1:2])
anova.cca(agua_rda, permutations = 10000)
RsquareAdj(agua_rda)

bacia_rda <- rda(X = assembleia_peixes_hell, Y =  bacia_PCs[,1:2])
anova.cca(bacia_rda, permutations = 10000)
RsquareAdj(bacia_rda)

urb_rda <- rda(X = assembleia_peixes_hell, Y =  urb)
anova.cca(urb_rda, permutations = 10000)
RsquareAdj(urb_rda)

#########################
only_estrutura_rda <- rda(X = assembleia_peixes_hell, Y =  estrutura_PCs[,1:2], Z = data.frame(agua_PCs[,1:2], bacia_PCs[,1:2], urb))
anova.cca(only_estrutura_rda, permutations = 10000)
RsquareAdj(only_estrutura_rda)

only_agua_rda <- rda(X = assembleia_peixes_hell, Y =  agua_PCs[,1:2], Z = data.frame(estrutura_PCs[,1:2], bacia_PCs[,1:2], urb))
anova.cca(only_agua_rda, permutations = 10000)
RsquareAdj(only_agua_rda)

only_bacia_rda <- rda(X = assembleia_peixes_hell, Y =  bacia_PCs[,1:2], Z = data.frame(agua_PCs[,1:2], estrutura_PCs[,1:2], urb))
anova.cca(only_bacia_rda, permutations = 10000)
RsquareAdj(only_bacia_rda)

only_urb_rda <- rda(X = assembleia_peixes_hell, Y =  urb, Z = data.frame(agua_PCs[,1:2], bacia_PCs[,1:2], estrutura_PCs[,1:2]))
anova.cca(only_urb_rda, permutations = 10000)
RsquareAdj(only_urb_rda)



pdf("plots/varpart_hellinger.pdf", height = 3.5, width = 5, pointsize = 5)
plot(varp2, Xnames = c(paste("Estrutura", round(RsquareAdj(estrutura_rda)$adj.r.squared,2)),
                       paste("Água",round(RsquareAdj(agua_rda)$adj.r.squared,2)),
                       paste("Bacia",round(RsquareAdj(bacia_rda)$adj.r.squared,2)),
                       paste("Urbanização",round(RsquareAdj(urb_rda)$adj.r.squared,2))))
dev.off()


