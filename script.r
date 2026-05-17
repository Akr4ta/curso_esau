# instalar pacotes caso não estejam instalados
lapply(c("bipartite", "econet", "car", "carData", "circlize", 
         "vegan", "igraph", "RColorBrewer", "betapart", "iNEXT", "devtools"),
  function(pkg) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }
)
install_github("paternogbc/ecodados")
# carregrar bibliotecas
lapply(c("bipartite", "econet", "car", "carData", "circlize", "vegan", "igraph", "RColorBrewer", "ecodados", "betapart", "iNEXT", "ggplot2"), library, character.only = TRUE)

# carregando tabelas
trophic_pon <- read.delim("../data/trophic_pon.txt", row.names=1)
trophic_bin <- read.delim("../data/trophic_bin.txt", row.names=1)
composicao_especies <- ecodados::composicao_anuros_div_taxonomica
precipitacao <- ecodados::precipitacao_div_taxonomica
composicao_especies

### Diversidade
# shannon
shannon_res <- vegan::diversity(composicao_especies, index = "shannon", MARGIN = 1)
shannon_res

# simpson
simpson_res <- vegan::diversity(composicao_especies, index = "simpson", MARGIN = 1) 
simpson_res

# pielou
Pielou <- shannon_res/log(specnumber(composicao_especies))
Pielou

# exemplo pratico
## Juntando todos os dados em um único data frame
dados_div <- data.frame(precipitacao$prec, shannon_res, simpson_res, Pielou)
## Renomeando as colunas
colnames(dados_div) <- c("Precipitacao", "Shannon", "Simpson", "Pielou")
dados_div

# shannon
anova_shan <- lm(Shannon ~ Precipitacao, data = dados_div)
anova(anova_shan)

# simpson
anova_simp <- lm(Simpson ~ Precipitacao, data = dados_div)
anova(anova_simp)

# pielou
anova_piel <- lm(Pielou ~ Precipitacao, data = dados_div)
anova(anova_piel)

### Dissimilaridade
## Binarios
# Transformando dados em presencia e ausência.
composicao_PA <- decostand(composicao_especies, method = "pa")
# Diversidade beta
resultado_PA <- beta.pair(composicao_PA, index.family = "sorensen")
resultado_PA$beta.sor

## Abundancia
# Diversidade beta para abundância
resultado_AB <- beta.pair.abund(composicao_especies, index.family = "bray")
resultado_AB$beta.bray

### Rarefação
## Dados
data("mite")
data("mite.xy")
coord <- mite.xy
colnames(coord) <- c("long", "lat")
data("mite.env")
agua <- mite.env[, 2]
dados_rarefacao <- ecodados::rarefacao_morcegos
rarefacao_repteis <- ecodados::rarefacao_repteis
rarefacao_anuros <- ecodados::rarefacao_anuros
dados_amostras <- ecodados::morcegos_rarefacao_amostras

## baseada no individuo
resultados_morcegos <- iNEXT(dados_rarefacao, q = 0, datatype = "abundance", endpoint = 800)

ggiNEXT(resultados_morcegos, type = 1) +
    geom_vline(xintercept = 166, lty = 2) +
    scale_linetype_discrete(labels = c("Interpolado", "Extrapolado")) +
    scale_colour_manual(values = c("darkorange", "darkorchid", "cyan4")) +
    scale_fill_manual(values = c("darkorange", "darkorchid", "cyan4")) +
    labs(x = "Número de indivíduos", y = " Riqueza de espécies") +
    tema_livro()

## baseado na amostras
lista_rarefacao <- list(Tenentes = dados_amostras[1:18, 1],
                        Talhadinho = dados_amostras[, 2],
                        Experimental = dados_amostras[1:16, 3])

# Análise
res_rarefacao_amostras <- iNEXT(lista_rarefacao, q = 0, 
                                datatype = "incidence_freq")

# Gráfico
ggiNEXT(res_rarefacao_amostras , type = 1) + 
    geom_vline(xintercept = 12, lty = 2) +
    scale_linetype_discrete(name = "Método", labels = c("Interpolado", "Extrapolado")) +
    scale_colour_manual(values = c("darkorange", "darkorchid", "cyan4")) +
    scale_fill_manual(values = c("darkorange", "darkorchid", "cyan4")) +
    labs(x = "Número de amostras", y = " Riqueza de espécies") +
    tema_livro()


### Redes de interacoes
NlResult <- networklevel(trophic_pon)
NlResult2 <- networklevel(trophic_bin)
NlResult
NlResult2
NlResult[c(1, 2, 4, 7, 8, 9, 11, 18, 19, 21, 30, 31, 42, 43, 44, 45, 49)]

## Graficos
# criando vetores
nome <- c("trofic_pon")
tabela <- c("trophic_pon")

# grafico bipartites
plotweb_deprecated(trophic_pon, method="normal", text.rot="90", labsize=0.9, col.low="deepskyblue4", col.high="gold1", col.interaction="cornsilk4")

# grafico circular
chordDiagram(as.matrix(trophic_pon), small.gap = 1,
            annotationTrack = trophic_pon, 
            preAllocateTracks = 1)

circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
sector.name = get.cell.meta.data("sector.index")

circos.text(CELL_META$xcenter, CELL_META$ylim[1],
            sector.name,
            facing = "clockwise",
            niceFacing = TRUE,
            cex = 0.42,
            adj = c(-0.1, 0.1))

circos.axis(h = "bottom",        
            major.at = NULL,    
            minor.ticks = 0,    
            labels = FALSE,     
            lwd = 1.5,           
            col = "black")      
}, bg.border = NA)








### outros calculos
# Dissimilaridade binaria X preciptacao
data.frame_PA <- data.frame(round(as.numeric(resultado_PA$beta.sor), 2),
                            round(as.numeric(resultado_PA$beta.sim), 2),
                            round(as.numeric(resultado_PA$beta.sne), 2))
colnames(data.frame_PA) <- c("Sorensen", "Simpson", "Aninhamento")
data.frame_PA
prec_dis <- vegdist(precipitacao, method = "euclidian")
dados_prec <- as.numeric(prec_dis) 
dados_dis <- data.frame(dados_prec, data.frame_PA)
dados_dis

anova_sore <-lm(Sorensen ~ dados_prec, data = dados_dis)
anova(anova_sore)

# Dissimilaridade binaria X preciptacao
## Data frame
# Vamos montar um data.frame com os resultados.
data.frame_AB <- data.frame(round(as.numeric(resultado_AB$beta.bray), 2),
                            round(as.numeric(resultado_AB$beta.bray.bal), 2),
                            round(as.numeric(resultado_AB$beta.bray.gra), 2))
colnames(data.frame_AB) <- c("Bray", "Balanceada", "Gradiente")
data.frame_AB
## Agora vamos juntar os resultados com a precipitação
dados_dis_AB <- data.frame(dados_prec, data.frame_AB)
# Avaliar a relação entre os valores de diversidade beta total e precipitação
anova_dis_AB <- lm(Bray ~ dados_prec, data = dados_dis_AB)
anova(anova_dis_AB)