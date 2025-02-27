Exemplo varpart_manyglm
================
Rodolfo Pelinson
2025-02-27

# Exemplo varpart_manyglm

``` r
dir<-("C:/Users/rodol/OneDrive/repos/Urban_fish_assemblages")
```

``` r
library(vegan)
library(mvabund)
```

Carredango planilhas

``` r
assembleia_peixes <- read.csv(paste(sep = "/",dir,"data/com_por_bacia.csv"), row.names = 1)
```

Removendo espécies com apenas duas ou menos presenças

``` r
source(paste(sep = "/",dir,"functions/remove_sp.R"))

ncol(assembleia_peixes)
```

    ## [1] 25

``` r
assembleia_peixes_rm <- remove_sp(assembleia_peixes, 2)
ncol(assembleia_peixes_rm)
```

    ## [1] 7

Redução absurda no número de espécies.

Agora carregando as planilhas ambientais:

``` r
agua_PCs <- read.csv(paste(sep = "/",dir,"data/pcas_amb/agua_PCs.csv"), row.names = 1)
estrutura_PCs <- read.csv(paste(sep = "/",dir,"data/pcas_amb/estrutura_PCs.csv"), row.names = 1)
bacia_PCs <- read.csv(paste(sep = "/",dir,"data/pcas_amb/bacia_PCs.csv"), row.names = 1)

delineamento <- read.csv(paste(sep = "/",dir,"data/delineamento.csv"))

ncol(agua_PCs)
```

    ## [1] 10

``` r
ncol(estrutura_PCs)
```

    ## [1] 29

``` r
ncol(bacia_PCs)
```

    ## [1] 7

Agora carregando planilha com distancias e gerando os filtros espaciais

``` r
library(ade4)
library(adespatial)
```

    ## Registered S3 methods overwritten by 'adegraphics':
    ##   method         from
    ##   biplot.dudi    ade4
    ##   kplot.foucart  ade4
    ##   kplot.mcoa     ade4
    ##   kplot.mfa      ade4
    ##   kplot.pta      ade4
    ##   kplot.sepan    ade4
    ##   kplot.statis   ade4
    ##   scatter.coa    ade4
    ##   scatter.dudi   ade4
    ##   scatter.nipals ade4
    ##   scatter.pco    ade4
    ##   score.acm      ade4
    ##   score.mix      ade4
    ##   score.pca      ade4
    ##   screeplot.dudi ade4

    ## Registered S3 method overwritten by 'spdep':
    ##   method   from
    ##   plot.mst ape

    ## Registered S3 methods overwritten by 'adespatial':
    ##   method             from       
    ##   plot.multispati    adegraphics
    ##   print.multispati   ade4       
    ##   summary.multispati ade4

    ## 
    ## Anexando pacote: 'adespatial'

    ## O seguinte objeto é mascarado por 'package:ade4':
    ## 
    ##     multispati

``` r
dist_euclid <- read.csv(paste(sep = "/",dir,"data/dist/Matriz_distancia_matriz_euclidiana.csv"), row.names = 1)
dist_euclid <- as.dist(dist_euclid)

dbmem_euclid <- dbmem(dist_euclid, thresh = NULL, MEM.autocor = c("positive", "non-null", "all", "negative"), store.listw = TRUE, silent = FALSE)
```

    ## Truncation level = 0.3268453 
    ## Time to compute dbMEMs = 0.010000  sec

``` r
ncol(dbmem_euclid)
```

    ## [1] 5

Agora executando a função forward selection para com “mvabund”. Essa
fução basicamente vai: 1 - Rodar um modelo da matriz resposta em função
de cada preditor da matriz preditora separadamente 2 - Reorganizar a
matriz preditora do preditor com maior R² para o menor 3 - Rodar o
modelo adicionando preditora por preditora em ordem do maior R² para o
menor, até que adicionar novos preditores não gere aumento no R² do
modelo, ou o modelo com o novo preditor não seja significativamente
melhor que o anterior.

``` r
source(paste(sep = "/",dir,"functions/R2_manyglm.R"))
source(paste(sep = "/",dir,"functions/forward_sel_manyglm.R"))
set.seed(1)
estrutura_FS <- forward_sel_manyglm(y = assembleia_peixes_rm, x = data.frame(estrutura_PCs[,1:5]), nBoot=4999) #considerarei apenas os 5 primeiros eixos
```

    ## Time elapsed: 0 hr 0 min 44 sec
    ## Time elapsed: 0 hr 1 min 9 sec
    ##     df.diff      Dev         R2      p
    ## PC2       1 18.12051 0.10244839 0.0442
    ## PC1       1 12.87992 0.08103233 0.2436

``` r
estrutura_FS
```

    ## $result
    ##     df.diff      Dev         R2      p
    ## PC2       1 18.12051 0.10244839 0.0442
    ## PC1       1 12.87992 0.08103233 0.2436
    ## 
    ## $new_x
    ##             PC2
    ## 1   0.212314269
    ## 2   0.053844806
    ## 3   0.065993448
    ## 4   0.049473900
    ## 5  -0.006413301
    ## 6   0.080531731
    ## 7  -0.038246342
    ## 8  -0.002264514
    ## 9   0.099227228
    ## 10  0.138894434
    ## 11 -0.495794982
    ## 12 -0.082350923
    ## 13  0.153498729
    ## 14  0.085186304
    ## 15  0.300209781
    ## 16 -0.117556474
    ## 17  0.120246741
    ## 18  0.046453084
    ## 19 -0.120779193
    ## 20  0.123092551
    ## 21  0.116386117
    ## 22 -0.252933677
    ## 23 -0.186495348
    ## 24  0.095774403
    ## 25  0.174598952
    ## 26  0.171467825
    ## 27 -0.254085256
    ## 28  0.111813496
    ## 29 -0.391091047
    ## 30 -0.250996741

``` r
agua_FS <- forward_sel_manyglm(y = assembleia_peixes_rm, x = data.frame(agua_PCs[,1:5]), nBoot=4999) #considerarei apenas os 5 primeiros eixos
```

    ## Time elapsed: 0 hr 0 min 43 sec
    ##     df.diff      Dev         R2      p
    ## PC1       1 14.68681 0.09402454 0.1238

``` r
agua_FS
```

    ## $result
    ##     df.diff      Dev         R2      p
    ## PC1       1 14.68681 0.09402454 0.1238
    ## 
    ## $new_x
    ## data frame com 0 coluna e 30 linhas

``` r
bacia_FS <- forward_sel_manyglm(y = assembleia_peixes_rm, x = data.frame(bacia_PCs[,1:5]), nBoot=4999) #considerarei apenas os 5 primeiros eixos
```

    ## Time elapsed: 0 hr 0 min 46 sec
    ## Time elapsed: 0 hr 1 min 1 sec
    ##     df.diff      Dev         R2      p
    ## PC1       1 19.25048 0.12376087 0.0382
    ## PC4       1 12.45575 0.08569951 0.2628

``` r
bacia_FS
```

    ## $result
    ##     df.diff      Dev         R2      p
    ## PC1       1 19.25048 0.12376087 0.0382
    ## PC4       1 12.45575 0.08569951 0.2628
    ## 
    ## $new_x
    ##              PC1
    ## 1  -2.129259e-01
    ## 2  -6.693936e-02
    ## 3  -3.346994e-02
    ## 4   1.608401e-01
    ## 5   2.586595e-02
    ## 6   5.583106e-02
    ## 7   2.617015e-02
    ## 8  -1.769776e-01
    ## 9   4.634636e-02
    ## 10  2.121586e-01
    ## 11  6.876881e-01
    ## 12 -2.751243e-01
    ## 13  8.173215e-02
    ## 14  6.380771e-02
    ## 15 -8.575594e-02
    ## 16 -1.452876e-01
    ## 17 -5.207655e-02
    ## 18 -8.771035e-02
    ## 19 -2.005602e-01
    ## 20  1.726985e-01
    ## 21  1.057264e-01
    ## 22 -1.608626e-01
    ## 23 -2.376238e-02
    ## 24 -6.194725e-07
    ## 25  1.123785e-02
    ## 26 -1.975821e-01
    ## 27 -1.586415e-01
    ## 28 -7.291164e-02
    ## 29  5.137353e-02
    ## 30  2.491122e-01

``` r
dbmem_FS <- forward_sel_manyglm(y = assembleia_peixes_rm, x = data.frame(dbmem_euclid), nBoot=4999) #considerarei apenas os 5 primeiros eixos
```

    ## Time elapsed: 0 hr 0 min 44 sec
    ## Time elapsed: 0 hr 1 min 1 sec
    ## Time elapsed: 0 hr 0 min 59 sec
    ##      df.diff      Dev         R2      p
    ## MEM1       1 32.57633 0.18166924 0.0022
    ## MEM4       1 27.28412 0.13556333 0.0344
    ## MEM5       1 12.43050 0.09125614 0.3024

``` r
dbmem_FS
```

    ## $result
    ##      df.diff      Dev         R2      p
    ## MEM1       1 32.57633 0.18166924 0.0022
    ## MEM4       1 27.28412 0.13556333 0.0344
    ## MEM5       1 12.43050 0.09125614 0.3024
    ## 
    ## $new_x
    ##              MEM1        MEM4
    ## b031  -1.20004534  1.01699034
    ## b034  -1.19853084  1.02269779
    ## b039  -1.05542434 -0.77027483
    ## b040  -1.05476486 -0.76402533
    ## b066  -0.85449392 -1.56406361
    ## b202   0.05283562  0.87818680
    ## b204   0.06175303  0.87137395
    ## b309  -0.85582492 -1.57695726
    ## b310  -1.19991376  1.00373333
    ## b320  -0.46020633 -0.40471479
    ## b321  -1.25820440  1.55143430
    ## b344  -1.24994474 -0.91317791
    ## b539  -1.05426850 -0.70120505
    ## b543  -0.20547094  0.36901961
    ## b545  -0.20086575  0.36306724
    ## b570   0.55480908 -0.78493513
    ## b574   0.34423238  0.25725574
    ## b578   0.35134770  0.27218578
    ## b579   1.04804256 -2.70736507
    ## b581  -0.44997943 -0.42275657
    ## b589   0.06257594  0.91910135
    ## b594   1.47239909 -1.23737215
    ## b620   1.83674142 -0.93406179
    ## b627   0.37160042  0.31423220
    ## b631   0.38108905  0.32190484
    ## b637   1.69869027  0.79841004
    ## b673   1.81910610  1.28035542
    ## b711   1.86486825  0.83166849
    ## ebbrc  0.10694226  0.06158481
    ## ebbsn  0.27090489  0.64770743

Aqui nenhuma variavel de estrutura ou agua foi significativa. Apenas o
primeiro eixo para bacia foi significativo, e 2 eixes espaciais foram
significativos.

# Realizando a partição de variância

Vou adotar 3 estratégias.

## Primeira estratégia

1 - Fazer a partição de variância usando o procedimento clássico, que é
com R² ajustado e usando as planilhas resultantes da forward selection.

Nesse caso, quando a forward selection não selecionou nenhuma variavel,
eu vou pegar só a primeira. O padrão é pegar todas, temos casos onde
temos mais preditores que comunidades.. então nem daria.

``` r
source(paste(sep = "/",dir,"functions/varpart_manyglm.R"))

#grupos de preditoras devem estar organizadas em listas (com nome)
predictors <- list(estrutura = data.frame(PC1 = estrutura_PCs[,1]),
                   bacia = bacia_FS$new_x,
                   agua = data.frame(PC1 = agua_PCs[,1]),
                   urbanizacao = data.frame(urbanizacao = delineamento$urbana),
                   MEMs = dbmem_FS$new_x)

varpart_peixes <-varpart_manyglm(resp = assembleia_peixes_rm, pred = predictors, DF_adj_r2 = TRUE)

#Valores pra comunidade toda
varpart_peixes$R2_fractions_com
```

    ##             R2_full_fraction R2_pure_fraction
    ## estrutura         0.09661211      0.021969521
    ## bacia             0.12376087     -0.009864169
    ## agua              0.09402454      0.006754420
    ## urbanizacao       0.13481463     -0.013799889
    ## MEMs              0.27209774      0.237765532

``` r
#Valores para cada espécie, valores cheios
varpart_peixes$R2_fractions_sp$R2_full_fraction
```

    ##                                  estrutura        bacia         agua
    ## Gymnotus_pantherinus          6.246974e-07 0.0221751332 5.301314e-02
    ## Phalloceros_harpagos          5.507945e-02 0.0567994530 1.011259e-01
    ## Phalloceros_reisi             7.302165e-03 0.0007081747 9.827428e-05
    ## Hollandichthys_multifasciatus 3.170013e-02 0.0048476545 1.807529e-03
    ## Astyanax_lacustris            2.250839e-01 0.2731285007 1.472663e-01
    ## Poecilia_reticulata           2.501511e-01 0.1427106962 1.844557e-01
    ## Hoplosternum_littorale        1.069674e-01 0.3659564746 1.704049e-01
    ##                               urbanizacao       MEMs
    ## Gymnotus_pantherinus          0.022710699 0.08796179
    ## Phalloceros_harpagos          0.059557698 0.39208652
    ## Phalloceros_reisi             0.013623428 0.25653306
    ## Hollandichthys_multifasciatus 0.009598844 0.01322681
    ## Astyanax_lacustris            0.243249704 0.12633640
    ## Poecilia_reticulata           0.319600985 0.77313132
    ## Hoplosternum_littorale        0.275361060 0.25540825

``` r
round(varpart_peixes$R2_fractions_sp$R2_full_fraction,4)
```

    ##                               estrutura  bacia   agua urbanizacao   MEMs
    ## Gymnotus_pantherinus             0.0000 0.0222 0.0530      0.0227 0.0880
    ## Phalloceros_harpagos             0.0551 0.0568 0.1011      0.0596 0.3921
    ## Phalloceros_reisi                0.0073 0.0007 0.0001      0.0136 0.2565
    ## Hollandichthys_multifasciatus    0.0317 0.0048 0.0018      0.0096 0.0132
    ## Astyanax_lacustris               0.2251 0.2731 0.1473      0.2432 0.1263
    ## Poecilia_reticulata              0.2502 0.1427 0.1845      0.3196 0.7731
    ## Hoplosternum_littorale           0.1070 0.3660 0.1704      0.2754 0.2554

``` r
#Valores para cada espécie, valores puros
varpart_peixes$R2_fractions_sp$R2_pure_fraction
```

    ##                                  estrutura        bacia         agua
    ## Gymnotus_pantherinus           0.102138676 -0.007670511  0.072154064
    ## Phalloceros_harpagos          -0.008746554  0.042594477  0.004288223
    ## Phalloceros_reisi              0.091934408 -0.012813017  0.001219894
    ## Hollandichthys_multifasciatus  0.054784912 -0.005205655  0.055873703
    ## Astyanax_lacustris            -0.028987185 -0.028987326 -0.028961463
    ## Poecilia_reticulata           -0.029275092 -0.029284750 -0.029275610
    ## Hoplosternum_littorale        -0.028062518 -0.027682404 -0.028017870
    ##                                urbanizacao       MEMs
    ## Gymnotus_pantherinus          -0.008651265 0.04450906
    ## Phalloceros_harpagos          -0.012805396 0.25565672
    ## Phalloceros_reisi             -0.010128948 0.27646278
    ## Hollandichthys_multifasciatus  0.001759931 0.06624216
    ## Astyanax_lacustris            -0.009667965 0.40326872
    ## Poecilia_reticulata           -0.029240768 0.35636486
    ## Hoplosternum_littorale        -0.027864812 0.26185441

``` r
round(varpart_peixes$R2_fractions_sp$R2_pure_fraction,4)
```

    ##                               estrutura   bacia    agua urbanizacao   MEMs
    ## Gymnotus_pantherinus             0.1021 -0.0077  0.0722     -0.0087 0.0445
    ## Phalloceros_harpagos            -0.0087  0.0426  0.0043     -0.0128 0.2557
    ## Phalloceros_reisi                0.0919 -0.0128  0.0012     -0.0101 0.2765
    ## Hollandichthys_multifasciatus    0.0548 -0.0052  0.0559      0.0018 0.0662
    ## Astyanax_lacustris              -0.0290 -0.0290 -0.0290     -0.0097 0.4033
    ## Poecilia_reticulata             -0.0293 -0.0293 -0.0293     -0.0292 0.3564
    ## Hoplosternum_littorale          -0.0281 -0.0277 -0.0280     -0.0279 0.2619

O Objeto resultante contém os R² para todas as combinações de grupos de
preditores, então é só usar esses valores para calcular frações
específicas.

Por exemplo, se quiser saber a fração compartilhada entre variaveis da
bacia e urbanização:

``` r
#Bacia sem urbanização
bacia_sem_urb <- varpart_peixes$R2_models$`bacia-urbanizacao` - varpart_peixes$R2_models$bacia
bacia_sem_urb
```

    ## [1] 0.05678596

``` r
#Urbanização sem bacia
urb_sem_bacia <- varpart_peixes$R2_models$`bacia-urbanizacao` - varpart_peixes$R2_models$urbanizacao
urb_sem_bacia
```

    ## [1] 0.0457322

``` r
#Fração compartilhada
varpart_peixes$R2_models$`bacia-urbanizacao` - (bacia_sem_urb + urb_sem_bacia)
```

    ## [1] 0.07802867

``` r
#Aqui o valor é superestimado porque as frações puras foram negativas.
```

Mesma coisa para cada espécie:

``` r
#Bacia sem urbanização
bacia_sem_urb <- varpart_peixes$R2_models_sp$bacia.urbanizacao - varpart_peixes$R2_models_sp$bacia
bacia_sem_urb
```

    ## [1] 0.024716397 0.121199350 0.015932154 0.004456576 0.016345631 0.173375391
    ## [7] 0.041476248

``` r
#Urbanização sem bacia
urb_sem_bacia <- varpart_peixes$R2_models_sp$bacia.urbanizacao - varpart_peixes$R2_models_sp$urbanizacao
urb_sem_bacia
```

    ## [1]  0.0241808317  0.1184411050  0.0030168999 -0.0002946141  0.0462244278
    ## [6] -0.0035148979  0.1320716634

``` r
#Fração compartilhada
frac_comp <- varpart_peixes$R2_models_sp$bacia.urbanizacao - (bacia_sem_urb + urb_sem_bacia)

#Adicionando nomes:
names(bacia_sem_urb)<- rownames(varpart_peixes$R2_models_sp)
names(urb_sem_bacia)<- rownames(varpart_peixes$R2_models_sp)
names(frac_comp)<- rownames(varpart_peixes$R2_models_sp)

bacia_sem_urb
```

    ##          Gymnotus_pantherinus          Phalloceros_harpagos 
    ##                   0.024716397                   0.121199350 
    ##             Phalloceros_reisi Hollandichthys_multifasciatus 
    ##                   0.015932154                   0.004456576 
    ##            Astyanax_lacustris           Poecilia_reticulata 
    ##                   0.016345631                   0.173375391 
    ##        Hoplosternum_littorale 
    ##                   0.041476248

``` r
urb_sem_bacia
```

    ##          Gymnotus_pantherinus          Phalloceros_harpagos 
    ##                  0.0241808317                  0.1184411050 
    ##             Phalloceros_reisi Hollandichthys_multifasciatus 
    ##                  0.0030168999                 -0.0002946141 
    ##            Astyanax_lacustris           Poecilia_reticulata 
    ##                  0.0462244278                 -0.0035148979 
    ##        Hoplosternum_littorale 
    ##                  0.1320716634

``` r
frac_comp
```

    ##          Gymnotus_pantherinus          Phalloceros_harpagos 
    ##                  -0.002005699                  -0.061641652 
    ##             Phalloceros_reisi Hollandichthys_multifasciatus 
    ##                  -0.002308725                   0.005142269 
    ##            Astyanax_lacustris           Poecilia_reticulata 
    ##                   0.226904073                   0.146225594 
    ##        Hoplosternum_littorale 
    ##                   0.233884811

Com exceção das frações compartilhadas, podemos calcular o valor de p
para qualquer fração! Por exemplo, vou calcular o valor de p para as
variaveis de bacia, a fração completa, a fração pura, e a descontando
apenas a urbanização.

Primeiro a completa:

``` r
anova.manyglm(varpart_peixes$model_null, varpart_peixes$models$bacia, nBoot=999, p.uni="adjusted")
```

    ## Time elapsed: 0 hr 0 min 9 sec

    ## Analysis of Deviance Table
    ## 
    ## varpart_peixes$model_null: resp_mv ~ 1
    ## varpart_peixes$models$bacia: resp_mv ~ PC1
    ## 
    ## Multivariate test:
    ##                             Res.Df Df.diff   Dev Pr(>Dev)  
    ## varpart_peixes$model_null       29                         
    ## varpart_peixes$models$bacia     28       1 19.25    0.038 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Univariate Tests:
    ##                             Gymnotus_pantherinus          Phalloceros_harpagos
    ##                                              Dev Pr(>Dev)                  Dev
    ## varpart_peixes$model_null                                                     
    ## varpart_peixes$models$bacia                0.495    0.775                1.248
    ##                                      Phalloceros_reisi         
    ##                             Pr(>Dev)               Dev Pr(>Dev)
    ## varpart_peixes$model_null                                      
    ## varpart_peixes$models$bacia    0.663             0.021    0.906
    ##                             Hollandichthys_multifasciatus         
    ##                                                       Dev Pr(>Dev)
    ## varpart_peixes$model_null                                         
    ## varpart_peixes$models$bacia                         0.111    0.906
    ##                             Astyanax_lacustris          Poecilia_reticulata
    ##                                            Dev Pr(>Dev)                 Dev
    ## varpart_peixes$model_null                                                  
    ## varpart_peixes$models$bacia              6.205    0.135               3.568
    ##                                      Hoplosternum_littorale         
    ##                             Pr(>Dev)                    Dev Pr(>Dev)
    ## varpart_peixes$model_null                                           
    ## varpart_peixes$models$bacia    0.273                  7.601    0.102
    ## Arguments:
    ##  Test statistics calculated assuming uncorrelated response (for faster computation) 
    ## P-value calculated using 999 iterations via PIT-trap resampling.

Agora a pura. Para a pura basicamente é só comparar o modelo com todos
os preditores, mas sem os de bacia, com o modelo com todos os
preditores.

``` r
anova.manyglm(varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`, varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`, nBoot=999, p.uni="adjusted")
```

    ## Time elapsed: 0 hr 0 min 12 sec

    ## Analysis of Deviance Table
    ## 
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`: resp_mv ~ PC1 + PC1.1 + urbanizacao + MEM1 + MEM4
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`: resp_mv ~ PC1 + PC1.1 + PC1.2 + urbanizacao + MEM1 + MEM4
    ## 
    ## Multivariate test:
    ##                                                               Res.Df Df.diff
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`           24        
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`     23       1
    ##                                                                 Dev Pr(>Dev)
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`                     
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs` 2.517     0.53
    ## 
    ## Univariate Tests:
    ##                                                               Gymnotus_pantherinus
    ##                                                                                Dev
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`                           
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`                0.031
    ##                                                                       
    ##                                                               Pr(>Dev)
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`               
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`    0.927
    ##                                                               Phalloceros_harpagos
    ##                                                                                Dev
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`                           
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`                2.368
    ##                                                                       
    ##                                                               Pr(>Dev)
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`               
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`    0.462
    ##                                                               Phalloceros_reisi
    ##                                                                             Dev
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`                        
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`             0.089
    ##                                                                       
    ##                                                               Pr(>Dev)
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`               
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`    0.927
    ##                                                               Hollandichthys_multifasciatus
    ##                                                                                         Dev
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`                                    
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`                         0.008
    ##                                                                       
    ##                                                               Pr(>Dev)
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`               
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`    0.927
    ##                                                               Astyanax_lacustris
    ##                                                                              Dev
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`                         
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`                  0
    ##                                                                       
    ##                                                               Pr(>Dev)
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`               
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`    0.927
    ##                                                               Poecilia_reticulata
    ##                                                                               Dev
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`                          
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`                   0
    ##                                                                       
    ##                                                               Pr(>Dev)
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`               
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`    0.927
    ##                                                               Hoplosternum_littorale
    ##                                                                                  Dev
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`                             
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`                   0.02
    ##                                                                       
    ##                                                               Pr(>Dev)
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`               
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`    0.927
    ## Arguments:
    ##  Test statistics calculated assuming uncorrelated response (for faster computation) 
    ## P-value calculated using 999 iterations via PIT-trap resampling.

Agora o valor de p das variaveis de bacia, descontando a urbanização.

``` r
anova.manyglm(varpart_peixes$models$`urbanizacao`, varpart_peixes$models$`bacia-urbanizacao`, nBoot=999, p.uni="adjusted")
```

    ## Time elapsed: 0 hr 0 min 11 sec

    ## Analysis of Deviance Table
    ## 
    ## varpart_peixes$models$urbanizacao: resp_mv ~ urbanizacao
    ## varpart_peixes$models$`bacia-urbanizacao`: resp_mv ~ PC1 + urbanizacao
    ## 
    ## Multivariate test:
    ##                                           Res.Df Df.diff  Dev Pr(>Dev)
    ## varpart_peixes$models$urbanizacao             28                      
    ## varpart_peixes$models$`bacia-urbanizacao`     27       1 8.74    0.419
    ## 
    ## Univariate Tests:
    ##                                           Gymnotus_pantherinus         
    ##                                                            Dev Pr(>Dev)
    ## varpart_peixes$models$urbanizacao                                      
    ## varpart_peixes$models$`bacia-urbanizacao`                0.588    0.899
    ##                                           Phalloceros_harpagos         
    ##                                                            Dev Pr(>Dev)
    ## varpart_peixes$models$urbanizacao                                      
    ## varpart_peixes$models$`bacia-urbanizacao`                2.946    0.493
    ##                                           Phalloceros_reisi         
    ##                                                         Dev Pr(>Dev)
    ## varpart_peixes$models$urbanizacao                                   
    ## varpart_peixes$models$`bacia-urbanizacao`              0.11    0.945
    ##                                           Hollandichthys_multifasciatus
    ##                                                                     Dev
    ## varpart_peixes$models$urbanizacao                                      
    ## varpart_peixes$models$`bacia-urbanizacao`                         0.001
    ##                                                    Astyanax_lacustris         
    ##                                           Pr(>Dev)                Dev Pr(>Dev)
    ## varpart_peixes$models$urbanizacao                                             
    ## varpart_peixes$models$`bacia-urbanizacao`    0.969              1.426    0.740
    ##                                           Poecilia_reticulata         
    ##                                                           Dev Pr(>Dev)
    ## varpart_peixes$models$urbanizacao                                     
    ## varpart_peixes$models$`bacia-urbanizacao`               0.246    0.945
    ##                                           Hoplosternum_littorale         
    ##                                                              Dev Pr(>Dev)
    ## varpart_peixes$models$urbanizacao                                        
    ## varpart_peixes$models$`bacia-urbanizacao`                  3.423    0.491
    ## Arguments:
    ##  Test statistics calculated assuming uncorrelated response (for faster computation) 
    ## P-value calculated using 999 iterations via PIT-trap resampling.

## Segunda estratégia

1 - Fazer a partição de variância usando um número padrão de variaveis
(2) para cada grupo de preditores, sem ajustar valores de R². A única
exceção é para urbanização que só tem 1 preditor. Vou considerar a
forward selection apenas para os dbMEMns, que selecionou apenas 2
preditores tb.

Eu prefiro assim, o ajuste dos valores de R² para adição de preditores
parece muito severo, fora as montes de valores negativos que aparecem,
fora o fato de que eu usei a mesma fórmula de ajuste do R² clássico, não
sei se ela se aplica a esse pseudo-R².

``` r
#grupos de preditoras devem estar organizadas em listas (com nome)
predictors <- list(estrutura = data.frame(PC1 = estrutura_PCs[,1:2]),
                   bacia = data.frame(PC1 = bacia_PCs[,1:2]),
                   agua = data.frame(PC1 = agua_PCs[,1:2]),
                   urbanizacao = data.frame(urbanizacao = delineamento$urbana),
                   MEMs = dbmem_FS$new_x)

varpart_peixes <-varpart_manyglm(resp = assembleia_peixes_rm, pred = predictors, DF_adj_r2 = FALSE)

#Valores pra comunidade toda
varpart_peixes$R2_fractions_com
```

    ##             R2_full_fraction R2_pure_fraction
    ## estrutura          0.1755151      0.088955307
    ## bacia              0.1585156      0.044049770
    ## agua               0.1635922      0.040224677
    ## urbanizacao        0.1348146      0.007270443
    ## MEMs               0.2818155      0.155494332

``` r
#Valores para cada espécie, valores cheios
varpart_peixes$R2_fractions_sp$R2_full_fraction
```

    ##                                estrutura      bacia        agua urbanizacao
    ## Gymnotus_pantherinus          0.13482064 0.02478269 0.122758095 0.022710699
    ## Phalloceros_harpagos          0.08167536 0.16668850 0.213237478 0.059557698
    ## Phalloceros_reisi             0.21000531 0.02172926 0.004440855 0.013623428
    ## Hollandichthys_multifasciatus 0.12013081 0.04426051 0.079138496 0.009598844
    ## Astyanax_lacustris            0.29322807 0.30673129 0.311204637 0.243249704
    ## Poecilia_reticulata           0.26628368 0.14613055 0.206198531 0.319600985
    ## Hoplosternum_littorale        0.12246151 0.39928665 0.208167305 0.275361060
    ##                                     MEMs
    ## Gymnotus_pantherinus          0.09110328
    ## Phalloceros_harpagos          0.40608961
    ## Phalloceros_reisi             0.26569496
    ## Hollandichthys_multifasciatus 0.01369919
    ## Astyanax_lacustris            0.13084842
    ## Poecilia_reticulata           0.80074315
    ## Hoplosternum_littorale        0.26452997

``` r
round(varpart_peixes$R2_fractions_sp$R2_full_fraction,4)
```

    ##                               estrutura  bacia   agua urbanizacao   MEMs
    ## Gymnotus_pantherinus             0.1348 0.0248 0.1228      0.0227 0.0911
    ## Phalloceros_harpagos             0.0817 0.1667 0.2132      0.0596 0.4061
    ## Phalloceros_reisi                0.2100 0.0217 0.0044      0.0136 0.2657
    ## Hollandichthys_multifasciatus    0.1201 0.0443 0.0791      0.0096 0.0137
    ## Astyanax_lacustris               0.2932 0.3067 0.3112      0.2432 0.1308
    ## Poecilia_reticulata              0.2663 0.1461 0.2062      0.3196 0.8007
    ## Hoplosternum_littorale           0.1225 0.3993 0.2082      0.2754 0.2645

``` r
#Valores para cada espécie, valores puros
varpart_peixes$R2_fractions_sp$R2_pure_fraction
```

    ##                                  estrutura        bacia         agua
    ## Gymnotus_pantherinus          1.297363e-01 9.695697e-02 1.031724e-01
    ## Phalloceros_harpagos          6.377210e-02 1.690581e-03 9.090093e-02
    ## Phalloceros_reisi             2.815826e-01 6.639394e-02 9.316207e-03
    ## Hollandichthys_multifasciatus 1.475256e-01 1.431858e-01 7.799416e-02
    ## Astyanax_lacustris            1.858871e-05 1.450787e-05 8.278103e-05
    ## Poecilia_reticulata           2.271945e-06 3.051030e-06 1.506970e-07
    ## Hoplosternum_littorale        4.966001e-05 1.035804e-04 1.060882e-04
    ##                                urbanizacao         MEMs
    ## Gymnotus_pantherinus          1.562423e-04 4.539012e-02
    ## Phalloceros_harpagos          4.715868e-04 5.386322e-01
    ## Phalloceros_reisi             2.674176e-02 3.333985e-01
    ## Hollandichthys_multifasciatus 2.350753e-02 1.708050e-01
    ## Astyanax_lacustris            5.651066e-07 6.082248e-05
    ## Poecilia_reticulata           1.374333e-05 1.458604e-04
    ## Hoplosternum_littorale        1.671014e-06 2.776444e-05

``` r
round(varpart_peixes$R2_fractions_sp$R2_pure_fraction,4)
```

    ##                               estrutura  bacia   agua urbanizacao   MEMs
    ## Gymnotus_pantherinus             0.1297 0.0970 0.1032      0.0002 0.0454
    ## Phalloceros_harpagos             0.0638 0.0017 0.0909      0.0005 0.5386
    ## Phalloceros_reisi                0.2816 0.0664 0.0093      0.0267 0.3334
    ## Hollandichthys_multifasciatus    0.1475 0.1432 0.0780      0.0235 0.1708
    ## Astyanax_lacustris               0.0000 0.0000 0.0001      0.0000 0.0001
    ## Poecilia_reticulata              0.0000 0.0000 0.0000      0.0000 0.0001
    ## Hoplosternum_littorale           0.0000 0.0001 0.0001      0.0000 0.0000

O Objeto resultante contém os R² para todas as combinações de grupos de
preditores, então é só usar esses valores para calcular frações
específicas.

Por exemplo, se quiser saber a fração compartilhada entre variaveis da
bacia e urbanização:

``` r
#Bacia sem urbanização
bacia_sem_urb <- varpart_peixes$R2_models$`bacia-urbanizacao` - varpart_peixes$R2_models$bacia
bacia_sem_urb
```

    ## [1] 0.135282

``` r
#Urbanização sem bacia
urb_sem_bacia <- varpart_peixes$R2_models$`bacia-urbanizacao` - varpart_peixes$R2_models$urbanizacao
urb_sem_bacia
```

    ## [1] 0.158983

``` r
#Fração compartilhada
varpart_peixes$R2_models$`bacia-urbanizacao` - (bacia_sem_urb + urb_sem_bacia) 
```

    ## [1] -0.0004673636

Mesma coisa para cada espécie:

``` r
#Bacia sem urbanização
bacia_sem_urb <- varpart_peixes$R2_models_sp$bacia.urbanizacao - varpart_peixes$R2_models_sp$bacia
bacia_sem_urb
```

    ## [1] 0.065725709 0.022342680 0.001182568 0.014345590 0.074410614 0.685766718
    ## [7] 0.083200085

``` r
#Urbanização sem bacia
urb_sem_bacia <- varpart_peixes$R2_models_sp$bacia.urbanizacao - varpart_peixes$R2_models_sp$urbanizacao 
urb_sem_bacia
```

    ## [1] 0.067797705 0.129473486 0.009288404 0.049007251 0.137892201 0.512296284
    ## [7] 0.207125679

``` r
#Fração compartilhada
frac_comp <- varpart_peixes$R2_models_sp$bacia.urbanizacao - (bacia_sem_urb + urb_sem_bacia)

#Adicionando nomes:
names(bacia_sem_urb)<- rownames(varpart_peixes$R2_models_sp)
names(urb_sem_bacia)<- rownames(varpart_peixes$R2_models_sp)
names(frac_comp)<- rownames(varpart_peixes$R2_models_sp)

bacia_sem_urb
```

    ##          Gymnotus_pantherinus          Phalloceros_harpagos 
    ##                   0.065725709                   0.022342680 
    ##             Phalloceros_reisi Hollandichthys_multifasciatus 
    ##                   0.001182568                   0.014345590 
    ##            Astyanax_lacustris           Poecilia_reticulata 
    ##                   0.074410614                   0.685766718 
    ##        Hoplosternum_littorale 
    ##                   0.083200085

``` r
urb_sem_bacia
```

    ##          Gymnotus_pantherinus          Phalloceros_harpagos 
    ##                   0.067797705                   0.129473486 
    ##             Phalloceros_reisi Hollandichthys_multifasciatus 
    ##                   0.009288404                   0.049007251 
    ##            Astyanax_lacustris           Poecilia_reticulata 
    ##                   0.137892201                   0.512296284 
    ##        Hoplosternum_littorale 
    ##                   0.207125679

``` r
frac_comp
```

    ##          Gymnotus_pantherinus          Phalloceros_harpagos 
    ##                  -0.043015011                   0.037215018 
    ##             Phalloceros_reisi Hollandichthys_multifasciatus 
    ##                   0.012440860                  -0.004746745 
    ##            Astyanax_lacustris           Poecilia_reticulata 
    ##                   0.168839090                  -0.366165733 
    ##        Hoplosternum_littorale 
    ##                   0.192160974

Com exceção das frações compartilhadas, podemos calcular o valor de p
para qualquer fração! Por exemplo, vou calcular o valor de p para as
variaveis de bacia, a fração completa, a fração pura, e a descontando
apenas a urbanização.

Primeiro a completa:

``` r
anova.manyglm(varpart_peixes$model_null, varpart_peixes$models$bacia, nBoot=999, p.uni="adjusted")
```

    ## Time elapsed: 0 hr 0 min 10 sec

    ## Analysis of Deviance Table
    ## 
    ## varpart_peixes$model_null: resp_mv ~ 1
    ## varpart_peixes$models$bacia: resp_mv ~ PC1.PC1 + PC1.PC2
    ## 
    ## Multivariate test:
    ##                             Res.Df Df.diff   Dev Pr(>Dev)
    ## varpart_peixes$model_null       29                       
    ## varpart_peixes$models$bacia     27       2 25.19    0.138
    ## 
    ## Univariate Tests:
    ##                             Gymnotus_pantherinus          Phalloceros_harpagos
    ##                                              Dev Pr(>Dev)                  Dev
    ## varpart_peixes$model_null                                                     
    ## varpart_peixes$models$bacia                0.554    0.880                3.821
    ##                                      Phalloceros_reisi         
    ##                             Pr(>Dev)               Dev Pr(>Dev)
    ## varpart_peixes$model_null                                      
    ## varpart_peixes$models$bacia    0.646             0.658    0.880
    ##                             Hollandichthys_multifasciatus         
    ##                                                       Dev Pr(>Dev)
    ## varpart_peixes$model_null                                         
    ## varpart_peixes$models$bacia                         1.031    0.880
    ##                             Astyanax_lacustris          Poecilia_reticulata
    ##                                            Dev Pr(>Dev)                 Dev
    ## varpart_peixes$model_null                                                  
    ## varpart_peixes$models$bacia              7.065    0.296               3.659
    ##                                      Hoplosternum_littorale         
    ##                             Pr(>Dev)                    Dev Pr(>Dev)
    ## varpart_peixes$model_null                                           
    ## varpart_peixes$models$bacia    0.646                    8.4    0.250
    ## Arguments:
    ##  Test statistics calculated assuming uncorrelated response (for faster computation) 
    ## P-value calculated using 999 iterations via PIT-trap resampling.

Agora a pura. Para a pura basicamente é só comparar o modelo com todos
os preditores, mas sem os de bacia, com o modelo com todos os
preditores.

``` r
anova.manyglm(varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`, varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`, nBoot=999, p.uni="adjusted")
```

    ## Time elapsed: 0 hr 0 min 13 sec

    ## Analysis of Deviance Table
    ## 
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`: resp_mv ~ PC1.PC1 + PC1.PC2 + PC1.PC1.1 + PC1.PC2.1 + urbanizacao + MEM1 + MEM4
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`: resp_mv ~ PC1.PC1 + PC1.PC2 + PC1.PC1.1 + PC1.PC2.1 + PC1.PC1.2 + PC1.PC2.2 + urbanizacao + MEM1 + MEM4
    ## 
    ## Multivariate test:
    ##                                                               Res.Df Df.diff
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`           22        
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`     20       2
    ##                                                                 Dev Pr(>Dev)
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`                     
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs` 12.11    0.375
    ## 
    ## Univariate Tests:
    ##                                                               Gymnotus_pantherinus
    ##                                                                                Dev
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`                           
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`                2.881
    ##                                                                       
    ##                                                               Pr(>Dev)
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`               
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`    0.566
    ##                                                               Phalloceros_harpagos
    ##                                                                                Dev
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`                           
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`                0.092
    ##                                                                       
    ##                                                               Pr(>Dev)
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`               
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`    0.566
    ##                                                               Phalloceros_reisi
    ##                                                                             Dev
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`                        
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`             4.857
    ##                                                                       
    ##                                                               Pr(>Dev)
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`               
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`    0.566
    ##                                                               Hollandichthys_multifasciatus
    ##                                                                                         Dev
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`                                    
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`                         4.276
    ##                                                                       
    ##                                                               Pr(>Dev)
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`               
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`    0.566
    ##                                                               Astyanax_lacustris
    ##                                                                              Dev
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`                         
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`              0.001
    ##                                                                       
    ##                                                               Pr(>Dev)
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`               
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`    0.566
    ##                                                               Poecilia_reticulata
    ##                                                                               Dev
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`                          
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`                   0
    ##                                                                       
    ##                                                               Pr(>Dev)
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`               
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`    0.566
    ##                                                               Hoplosternum_littorale
    ##                                                                                  Dev
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`                             
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`                  0.004
    ##                                                                       
    ##                                                               Pr(>Dev)
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`               
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`    0.566
    ## Arguments:
    ##  Test statistics calculated assuming uncorrelated response (for faster computation) 
    ## P-value calculated using 999 iterations via PIT-trap resampling.

Agora o valor de p das variaveis de bacia, descontando a urbanização.

``` r
anova.manyglm(varpart_peixes$models$`urbanizacao`, varpart_peixes$models$`bacia-urbanizacao`, nBoot=999, p.uni="adjusted")
```

    ## Time elapsed: 0 hr 0 min 12 sec

    ## Analysis of Deviance Table
    ## 
    ## varpart_peixes$models$urbanizacao: resp_mv ~ urbanizacao
    ## varpart_peixes$models$`bacia-urbanizacao`: resp_mv ~ PC1.PC1 + PC1.PC2 + urbanizacao
    ## 
    ## Multivariate test:
    ##                                           Res.Df Df.diff   Dev Pr(>Dev)  
    ## varpart_peixes$models$urbanizacao             28                         
    ## varpart_peixes$models$`bacia-urbanizacao`     26       2 37.76    0.067 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Univariate Tests:
    ##                                           Gymnotus_pantherinus         
    ##                                                            Dev Pr(>Dev)
    ## varpart_peixes$models$urbanizacao                                      
    ## varpart_peixes$models$`bacia-urbanizacao`                1.567    0.832
    ##                                           Phalloceros_harpagos         
    ##                                                            Dev Pr(>Dev)
    ## varpart_peixes$models$urbanizacao                                      
    ## varpart_peixes$models$`bacia-urbanizacao`                3.062    0.735
    ##                                           Phalloceros_reisi         
    ##                                                         Dev Pr(>Dev)
    ## varpart_peixes$models$urbanizacao                                   
    ## varpart_peixes$models$`bacia-urbanizacao`             0.284    0.864
    ##                                           Hollandichthys_multifasciatus
    ##                                                                     Dev
    ## varpart_peixes$models$urbanizacao                                      
    ## varpart_peixes$models$`bacia-urbanizacao`                         1.152
    ##                                                    Astyanax_lacustris         
    ##                                           Pr(>Dev)                Dev Pr(>Dev)
    ## varpart_peixes$models$urbanizacao                                             
    ## varpart_peixes$models$`bacia-urbanizacao`    0.832              3.603    0.735
    ##                                           Poecilia_reticulata         
    ##                                                           Dev Pr(>Dev)
    ## varpart_peixes$models$urbanizacao                                     
    ## varpart_peixes$models$`bacia-urbanizacao`              23.136    0.021
    ##                                           Hoplosternum_littorale         
    ##                                                              Dev Pr(>Dev)
    ## varpart_peixes$models$urbanizacao                                        
    ## varpart_peixes$models$`bacia-urbanizacao`                  4.958    0.718
    ## Arguments:
    ##  Test statistics calculated assuming uncorrelated response (for faster computation) 
    ## P-value calculated using 999 iterations via PIT-trap resampling.

## Terceira estratégia

1 - Fazer a partição de variância usando um número de PCs que indiquem
pelo menos 50% da variação nos dados para cada grupo de preditores, sem
ajustar valores de R². A única exceção é para urbanização que só tem 1
preditor. Vou considerar a forward selection apenas para os dbMEMns, que
selecionou apenas 2 preditores tb.

O ajuste dos valores de R² para adição de preditores parece muito
severo, fora as montes de valores negativos que aparecem, fora o fato de
que eu usei a mesma fórmula de ajuste do R² clássico, não sei se ela se
aplica a esse pseudo-R².

``` r
autovar_estrutura <- read.csv(paste(sep = "/",dir,"data/pcas_amb/estrutura_autovalores.csv"), row.names = 1)
autovar_agua <- read.csv(paste(sep = "/",dir,"data/pcas_amb/agua_autovalores.csv"), row.names = 1)
autovar_bacia <- read.csv(paste(sep = "/",dir,"data/pcas_amb/bacia_autovalores.csv"), row.names = 1)

autovar_estrutura$importance[1]
```

    ## [1] 0.17

``` r
autovar_estrutura$importance[1] + autovar_estrutura$importance[2]
```

    ## [1] 0.29

``` r
autovar_estrutura$importance[1] + autovar_estrutura$importance[2] + autovar_estrutura$importance[3]
```

    ## [1] 0.38

``` r
autovar_estrutura$importance[1] + autovar_estrutura$importance[2] + autovar_estrutura$importance[3] + autovar_estrutura$importance[4] + autovar_estrutura$importance[5]
```

    ## [1] 0.52

``` r
autovar_agua$importance[1]
```

    ## [1] 0.74

``` r
autovar_bacia$importance[1] + autovar_bacia$importance[2]
```

    ## [1] 0.63

``` r
#grupos de preditoras devem estar organizadas em listas (com nome)
predictors <- list(estrutura = data.frame(PC1 = estrutura_PCs[,1:5]),
                   bacia = data.frame(PC1 = bacia_PCs[,1:2]),
                   agua = data.frame(PC1 = agua_PCs[,1]),
                   urbanizacao = data.frame(urbanizacao = delineamento$urbana),
                   MEMs = dbmem_FS$new_x)

varpart_peixes <-varpart_manyglm(resp = assembleia_peixes_rm, pred = predictors, DF_adj_r2 = FALSE)

#Valores pra comunidade toda
varpart_peixes$R2_fractions_com
```

    ##             R2_full_fraction R2_pure_fraction
    ## estrutura         0.26246384       0.24820834
    ## bacia             0.15851564       0.08721191
    ## agua              0.09402454       0.15293137
    ## urbanizacao       0.13481463       0.01614095
    ## MEMs              0.28181551       0.21155125

``` r
#Valores para cada espécie, valores cheios
varpart_peixes$R2_fractions_sp$R2_full_fraction
```

    ##                               estrutura      bacia         agua urbanizacao
    ## Gymnotus_pantherinus          0.2634084 0.02478269 5.301314e-02 0.022710699
    ## Phalloceros_harpagos          0.1456027 0.16668850 1.011259e-01 0.059557698
    ## Phalloceros_reisi             0.2713583 0.02172926 9.827428e-05 0.013623428
    ## Hollandichthys_multifasciatus 0.1493707 0.04426051 1.807529e-03 0.009598844
    ## Astyanax_lacustris            0.4968391 0.30673129 1.472663e-01 0.243249704
    ## Poecilia_reticulata           0.3116659 0.14613055 1.844557e-01 0.319600985
    ## Hoplosternum_littorale        0.1990019 0.39928665 1.704049e-01 0.275361060
    ##                                     MEMs
    ## Gymnotus_pantherinus          0.09110328
    ## Phalloceros_harpagos          0.40608961
    ## Phalloceros_reisi             0.26569496
    ## Hollandichthys_multifasciatus 0.01369919
    ## Astyanax_lacustris            0.13084842
    ## Poecilia_reticulata           0.80074315
    ## Hoplosternum_littorale        0.26452997

``` r
round(varpart_peixes$R2_fractions_sp$R2_full_fraction,4)
```

    ##                               estrutura  bacia   agua urbanizacao   MEMs
    ## Gymnotus_pantherinus             0.2634 0.0248 0.0530      0.0227 0.0911
    ## Phalloceros_harpagos             0.1456 0.1667 0.1011      0.0596 0.4061
    ## Phalloceros_reisi                0.2714 0.0217 0.0001      0.0136 0.2657
    ## Hollandichthys_multifasciatus    0.1494 0.0443 0.0018      0.0096 0.0137
    ## Astyanax_lacustris               0.4968 0.3067 0.1473      0.2432 0.1308
    ## Poecilia_reticulata              0.3117 0.1461 0.1845      0.3196 0.8007
    ## Hoplosternum_littorale           0.1990 0.3993 0.1704      0.2754 0.2645

``` r
#Valores para cada espécie, valores puros
varpart_peixes$R2_fractions_sp$R2_pure_fraction
```

    ##                                  estrutura        bacia         agua
    ## Gymnotus_pantherinus          6.266611e-01 2.123536e-04 4.886121e-01
    ## Phalloceros_harpagos          6.289716e-02 5.186783e-04 1.120441e-05
    ## Phalloceros_reisi             3.896903e-01 6.867843e-02 6.591869e-02
    ## Hollandichthys_multifasciatus 6.580378e-01 5.410221e-01 5.159639e-01
    ## Astyanax_lacustris            4.901822e-05 8.024371e-06 2.179142e-06
    ## Poecilia_reticulata           3.622319e-06 5.795967e-06 2.411793e-07
    ## Hoplosternum_littorale        1.192691e-04 3.791940e-05 1.117607e-05
    ##                                urbanizacao         MEMs
    ## Gymnotus_pantherinus          8.798271e-05 2.063159e-03
    ## Phalloceros_harpagos          4.188511e-06 5.587719e-01
    ## Phalloceros_reisi             6.589317e-02 2.773461e-01
    ## Hollandichthys_multifasciatus 4.693446e-02 6.425748e-01
    ## Astyanax_lacustris            4.104759e-05 4.911133e-05
    ## Poecilia_reticulata           1.348217e-05 4.648027e-05
    ## Hoplosternum_littorale        1.230039e-05 7.133657e-06

``` r
round(varpart_peixes$R2_fractions_sp$R2_pure_fraction,4)
```

    ##                               estrutura  bacia   agua urbanizacao   MEMs
    ## Gymnotus_pantherinus             0.6267 0.0002 0.4886      0.0001 0.0021
    ## Phalloceros_harpagos             0.0629 0.0005 0.0000      0.0000 0.5588
    ## Phalloceros_reisi                0.3897 0.0687 0.0659      0.0659 0.2773
    ## Hollandichthys_multifasciatus    0.6580 0.5410 0.5160      0.0469 0.6426
    ## Astyanax_lacustris               0.0000 0.0000 0.0000      0.0000 0.0000
    ## Poecilia_reticulata              0.0000 0.0000 0.0000      0.0000 0.0000
    ## Hoplosternum_littorale           0.0001 0.0000 0.0000      0.0000 0.0000

O Objeto resultante contém os R² para todas as combinações de grupos de
preditores, então é só usar esses valores para calcular frações
específicas.

Por exemplo, se quiser saber a fração compartilhada entre variaveis da
bacia e urbanização:

``` r
#Bacia sem urbanização
bacia_sem_urb <- varpart_peixes$R2_models$`bacia-urbanizacao` - varpart_peixes$R2_models$bacia
bacia_sem_urb
```

    ## [1] 0.135282

``` r
#Urbanização sem bacia
urb_sem_bacia <- varpart_peixes$R2_models$`bacia-urbanizacao` - varpart_peixes$R2_models$urbanizacao
urb_sem_bacia
```

    ## [1] 0.158983

``` r
#Fração compartilhada
varpart_peixes$R2_models$`bacia-urbanizacao` - (bacia_sem_urb + urb_sem_bacia)
```

    ## [1] -0.0004673636

Mesma coisa para cada espécie:

``` r
#Bacia sem urbanização
bacia_sem_urb <- varpart_peixes$R2_models_sp$bacia.urbanizacao - varpart_peixes$R2_models_sp$bacia
bacia_sem_urb
```

    ## [1] 0.065725709 0.022342680 0.001182568 0.014345590 0.074410614 0.685766718
    ## [7] 0.083200085

``` r
#Urbanização sem bacia
urb_sem_bacia <- varpart_peixes$R2_models_sp$bacia.urbanizacao - varpart_peixes$R2_models_sp$urbanizacao 
urb_sem_bacia
```

    ## [1] 0.067797705 0.129473486 0.009288404 0.049007251 0.137892201 0.512296284
    ## [7] 0.207125679

``` r
#Fração compartilhada
frac_comp <- varpart_peixes$R2_models_sp$bacia.urbanizacao - (bacia_sem_urb + urb_sem_bacia)

#Adicionando nomes:
names(bacia_sem_urb)<- rownames(varpart_peixes$R2_models_sp)
names(urb_sem_bacia)<- rownames(varpart_peixes$R2_models_sp)
names(frac_comp)<- rownames(varpart_peixes$R2_models_sp)

bacia_sem_urb
```

    ##          Gymnotus_pantherinus          Phalloceros_harpagos 
    ##                   0.065725709                   0.022342680 
    ##             Phalloceros_reisi Hollandichthys_multifasciatus 
    ##                   0.001182568                   0.014345590 
    ##            Astyanax_lacustris           Poecilia_reticulata 
    ##                   0.074410614                   0.685766718 
    ##        Hoplosternum_littorale 
    ##                   0.083200085

``` r
urb_sem_bacia
```

    ##          Gymnotus_pantherinus          Phalloceros_harpagos 
    ##                   0.067797705                   0.129473486 
    ##             Phalloceros_reisi Hollandichthys_multifasciatus 
    ##                   0.009288404                   0.049007251 
    ##            Astyanax_lacustris           Poecilia_reticulata 
    ##                   0.137892201                   0.512296284 
    ##        Hoplosternum_littorale 
    ##                   0.207125679

``` r
frac_comp
```

    ##          Gymnotus_pantherinus          Phalloceros_harpagos 
    ##                  -0.043015011                   0.037215018 
    ##             Phalloceros_reisi Hollandichthys_multifasciatus 
    ##                   0.012440860                  -0.004746745 
    ##            Astyanax_lacustris           Poecilia_reticulata 
    ##                   0.168839090                  -0.366165733 
    ##        Hoplosternum_littorale 
    ##                   0.192160974

Com exceção das frações compartilhadas, podemos calcular o valor de p
para qualquer fração! Por exemplo, vou calcular o valor de p para as
variaveis de bacia, a fração completa, a fração pura, e a descontando
apenas a urbanização.

Primeiro a completa:

``` r
anova.manyglm(varpart_peixes$model_null, varpart_peixes$models$bacia, nBoot=999, p.uni="adjusted")
```

    ## Time elapsed: 0 hr 0 min 10 sec

    ## Analysis of Deviance Table
    ## 
    ## varpart_peixes$model_null: resp_mv ~ 1
    ## varpart_peixes$models$bacia: resp_mv ~ PC1.PC1 + PC1.PC2
    ## 
    ## Multivariate test:
    ##                             Res.Df Df.diff   Dev Pr(>Dev)
    ## varpart_peixes$model_null       29                       
    ## varpart_peixes$models$bacia     27       2 25.19    0.144
    ## 
    ## Univariate Tests:
    ##                             Gymnotus_pantherinus          Phalloceros_harpagos
    ##                                              Dev Pr(>Dev)                  Dev
    ## varpart_peixes$model_null                                                     
    ## varpart_peixes$models$bacia                0.554    0.903                3.821
    ##                                      Phalloceros_reisi         
    ##                             Pr(>Dev)               Dev Pr(>Dev)
    ## varpart_peixes$model_null                                      
    ## varpart_peixes$models$bacia    0.658             0.658    0.903
    ##                             Hollandichthys_multifasciatus         
    ##                                                       Dev Pr(>Dev)
    ## varpart_peixes$model_null                                         
    ## varpart_peixes$models$bacia                         1.031    0.903
    ##                             Astyanax_lacustris          Poecilia_reticulata
    ##                                            Dev Pr(>Dev)                 Dev
    ## varpart_peixes$model_null                                                  
    ## varpart_peixes$models$bacia              7.065    0.313               3.659
    ##                                      Hoplosternum_littorale         
    ##                             Pr(>Dev)                    Dev Pr(>Dev)
    ## varpart_peixes$model_null                                           
    ## varpart_peixes$models$bacia    0.658                    8.4    0.275
    ## Arguments:
    ##  Test statistics calculated assuming uncorrelated response (for faster computation) 
    ## P-value calculated using 999 iterations via PIT-trap resampling.

Agora a pura. Para a pura basicamente é só comparar o modelo com todos
os preditores, mas sem os de bacia, com o modelo com todos os
preditores.

``` r
anova.manyglm(varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`, varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`, nBoot=999, p.uni="adjusted")
```

    ## Time elapsed: 0 hr 0 min 12 sec

    ## Analysis of Deviance Table
    ## 
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`: resp_mv ~ PC1.PC1 + PC1.PC2 + PC1.PC3 + PC1.PC4 + PC1.PC5 + PC1 + urbanizacao + MEM1 + MEM4
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`: resp_mv ~ PC1.PC1 + PC1.PC2 + PC1.PC3 + PC1.PC4 + PC1.PC5 + PC1.PC1.1 + PC1.PC2.1 + PC1 + urbanizacao + MEM1 + MEM4
    ## 
    ## Multivariate test:
    ##                                                               Res.Df Df.diff
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`           20        
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`     18       2
    ##                                                                 Dev Pr(>Dev)
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`                     
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs` 28.75    0.197
    ## 
    ## Univariate Tests:
    ##                                                               Gymnotus_pantherinus
    ##                                                                                Dev
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`                           
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`                0.012
    ##                                                                       
    ##                                                               Pr(>Dev)
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`               
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`    0.483
    ##                                                               Phalloceros_harpagos
    ##                                                                                Dev
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`                           
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`                0.028
    ##                                                                       
    ##                                                               Pr(>Dev)
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`               
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`    0.483
    ##                                                               Phalloceros_reisi
    ##                                                                             Dev
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`                        
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`             5.843
    ##                                                                       
    ##                                                               Pr(>Dev)
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`               
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`    0.462
    ##                                                               Hollandichthys_multifasciatus
    ##                                                                                         Dev
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`                                    
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`                        22.864
    ##                                                                       
    ##                                                               Pr(>Dev)
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`               
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`    0.206
    ##                                                               Astyanax_lacustris
    ##                                                                              Dev
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`                         
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`                  0
    ##                                                                       
    ##                                                               Pr(>Dev)
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`               
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`    0.605
    ##                                                               Poecilia_reticulata
    ##                                                                               Dev
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`                          
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`                   0
    ##                                                                       
    ##                                                               Pr(>Dev)
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`               
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`    0.605
    ##                                                               Hoplosternum_littorale
    ##                                                                                  Dev
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`                             
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`                  0.001
    ##                                                                       
    ##                                                               Pr(>Dev)
    ## varpart_peixes$models$`estrutura-agua-urbanizacao-MEMs`               
    ## varpart_peixes$models$`estrutura-bacia-agua-urbanizacao-MEMs`    0.483
    ## Arguments:
    ##  Test statistics calculated assuming uncorrelated response (for faster computation) 
    ## P-value calculated using 999 iterations via PIT-trap resampling.

Agora o valor de p das variaveis de bacia, descontando a urbanização.

``` r
anova.manyglm(varpart_peixes$models$`urbanizacao`, varpart_peixes$models$`bacia-urbanizacao`, nBoot=999, p.uni="adjusted")
```

    ## Time elapsed: 0 hr 0 min 12 sec

    ## Analysis of Deviance Table
    ## 
    ## varpart_peixes$models$urbanizacao: resp_mv ~ urbanizacao
    ## varpart_peixes$models$`bacia-urbanizacao`: resp_mv ~ PC1.PC1 + PC1.PC2 + urbanizacao
    ## 
    ## Multivariate test:
    ##                                           Res.Df Df.diff   Dev Pr(>Dev)  
    ## varpart_peixes$models$urbanizacao             28                         
    ## varpart_peixes$models$`bacia-urbanizacao`     26       2 37.76     0.07 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Univariate Tests:
    ##                                           Gymnotus_pantherinus         
    ##                                                            Dev Pr(>Dev)
    ## varpart_peixes$models$urbanizacao                                      
    ## varpart_peixes$models$`bacia-urbanizacao`                1.567    0.830
    ##                                           Phalloceros_harpagos         
    ##                                                            Dev Pr(>Dev)
    ## varpart_peixes$models$urbanizacao                                      
    ## varpart_peixes$models$`bacia-urbanizacao`                3.062    0.706
    ##                                           Phalloceros_reisi         
    ##                                                         Dev Pr(>Dev)
    ## varpart_peixes$models$urbanizacao                                   
    ## varpart_peixes$models$`bacia-urbanizacao`             0.284    0.830
    ##                                           Hollandichthys_multifasciatus
    ##                                                                     Dev
    ## varpart_peixes$models$urbanizacao                                      
    ## varpart_peixes$models$`bacia-urbanizacao`                         1.152
    ##                                                    Astyanax_lacustris         
    ##                                           Pr(>Dev)                Dev Pr(>Dev)
    ## varpart_peixes$models$urbanizacao                                             
    ## varpart_peixes$models$`bacia-urbanizacao`    0.830              3.603    0.706
    ##                                           Poecilia_reticulata         
    ##                                                           Dev Pr(>Dev)
    ## varpart_peixes$models$urbanizacao                                     
    ## varpart_peixes$models$`bacia-urbanizacao`              23.136    0.022
    ##                                           Hoplosternum_littorale         
    ##                                                              Dev Pr(>Dev)
    ## varpart_peixes$models$urbanizacao                                        
    ## varpart_peixes$models$`bacia-urbanizacao`                  4.958    0.644
    ## Arguments:
    ##  Test statistics calculated assuming uncorrelated response (for faster computation) 
    ## P-value calculated using 999 iterations via PIT-trap resampling.
