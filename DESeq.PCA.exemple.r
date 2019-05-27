##Exemple d'analyse pour les RNAseq
#Référence: https://www.huber.embl.de/users/klaus/Teaching/DESeq2Predoc2014.html

#Installation des librairies
# faire une function d'installation qui vérifie la présence d'un package, installe au besoin, puis load la library
#Référence: https://gist.github.com/smithdanielle/9913897 avec modification
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, library, character.only = TRUE)
}
#Établir la liste de package à installer
pkg.rnaseq <- c("BiocStyle", "rmarkdown", "geneplotter", "ggplot2", "plyr",
                "LSD", "DESeq2", "gplots", "RColorBrewer", "stringr", "topGO",
                "genefilter", "biomaRt", "dplyr", "EDASeq", "fdrtool", "org.Mm.eg.db")
#Utiliser la fonction avec la liste de package
check.packages(pkg = pkg.rnaseq)

# problème d'installaion avec les packages: 
#'BiocStyle’, ‘geneplotter’, ‘DESeq2’, ‘topGO’, ‘genefilter’, ‘biomaRt’, ‘EDASeq’, ‘org.Mm.eg.db’
# Teste l'installation avec Bioconducteur
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
#Installation réussi sauf pour un update des pakages: "cluster" et "nlme"
BiocManager::install(c("BiocStyle", "geneplotter",  "topGO", "genefilter", "biomaRt", "EDASeq", "org.Mm.eg.db"))
library(DESeq2)
library(BiocStyle)
library(geneplotter)
library(topGO)
library(genefilter)
library(biomaRt)
library(EDASeq)
library(org.Mm.eg.db)
###Possibilité de réutiliser le code de la function check.package en le modifiant légèrement
#set working directory = done in opening the project in version control

#Load les données tests
load(url("http://www-huber.embl.de/users/klaus/geneCounts.RData"))
#Regarde l'objet DESeq2Table
#Comment créer l'objet? À voir... 
DESeq2Table
head(assay(DESeq2Table))
colSums(assay(DESeq2Table))
colData(DESeq2Table)
mcols(rowData(DESeq2Table))

### Correction of an annotation error, data are inverted
con <- as.character(colData(DESeq2Table)$condition)
con[2] <- "Del_8_17_homo"
con[7] <- "WT"
colData(DESeq2Table)$condition <- factor(con)

##??Définition d'un design exp et d'une formule qui est informative sur les groupes/éch. 
##Et qui est associé avec l'objet DESeqDataSet

##Quality control et Normalization 
#Compte le nombre de gènes qui sont différents de zéro
GeneCounts <- counts(DESeq2Table)
idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)})
sum(idx.nz)

##Pour la suite, on utilise un échantillon random de la table. 
#renomme la table avec les valeurs non-zéro 
nz.counts <- subset(GeneCounts, idx.nz)
#prend un échamtillon (sample)
sam <- sample(dim(nz.counts)[1], 5)
nz.counts[sam, ]

#Normalisation des données 
#1ere condition = renomme les déletions 'Del' et assure que la comparaison est toujours WT - Del
colData(DESeq2Table)$condition <- factor(colData(DESeq2Table)$condition, 
                                         levels = c("WT", "Del"))
colData(DESeq2Table)$condition <- factor(ifelse(is.na(colData(DESeq2Table)$condition),  
                                                "Del", "WT"), levels = c("WT", "Del"))
#Estimation of size factor
DESeq2Table <- estimateSizeFactors(DESeq2Table)
sizeFactors(DESeq2Table)
#Visualisation des densitées de comptes pour chaque échantillons
##Ne fonctionne pas, à revoir plus tard 
multidensity(counts(DESeq2Table, normalized = T)[idx.nz ,],
              xlab="mean counts", xlim=c(0, 1000))
#2e visualisation 
##Ne marche pas nonplus, mais mon objet DESeq n'es pas ID comme il faut je crois
multiecdf( counts(DESeq2Table, normalized = T)[idx.nz ,],
           xlab="mean counts", xlim=c(0, 1000))

#Autre visualisation, celle-ci fonctionnent bien
pdf("pairwiseMAs.pdf")
MA.idx = t(combn(1:8, 2))
for( i in  seq_along( MA.idx[,1])){ 
  MDPlot(counts(DESeq2Table, normalized = T)[idx.nz ,], 
         c(MA.idx[i,1],MA.idx[i,2]), 
         main = paste( colnames(DESeq2Table)[MA.idx[i,1]], " vs ",
                       colnames(DESeq2Table)[MA.idx[i,2]] ), ylim = c(-3,3))
}
dev.off()


###PCA and Heatmap, yeah!! 
#Normalisation logaritmique régularisé des données 
#Permet d'éviter le biais des petits ET très grand comptes
rld <- rlogTransformation(DESeq2Table, blind=TRUE)
#Calcule la distance entre le échantillons (transpose la matrice pour que la fonction marche)
distsRL <- dist(t(assay(rld)))
#Transforme les distance en matrice
mat <- as.matrix(distsRL)
rownames(mat) <-  colData(rld)$condition
colnames(mat) <-  colData(rld)$sampleNO
#Plot le heatmap
hmcol <- colorRampPalette(brewer.pal(9, "Blues"))(255)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))

#Plot un PCA des 500 top gènes
ntop = 500
Pvars <- rowVars(assay(rld))
#Sélection des 500 gènes avec la plus haute variance
select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, 
                                                      length(Pvars)))]
PCA <- prcomp(t(assay(rld)[select, ]), scale = F)
#Arrondi le % de variation
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
#Plot le graph avec ggplot
#Création de l'objet de graphqiue
dataGG = data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], 
                    PC3 = PCA$x[,3], PC4 = PCA$x[,4], 
                    sampleNO = colData(rld)$sampleNO,
                    condition = colData(rld)$condition)
#Plot le graph
(qplot(PC1, PC2, data = dataGG, color =  condition, 
       main = "PC1 vs PC2, top variable genes", size = I(6))
  + labs(x = paste0("PC1, VarExp:", round(percentVar[1],4)),
         y = paste0("PC2, VarExp:", round(percentVar[2],4)))
  + scale_colour_brewer(type="qual", palette=2)
)

###Retraits d'outliers au besoin 
#Se fait selon les cas et la valeurs choisi est variable et arbittraire
outliers <- as.character(subset(colnames(DESeq2Table), dataGG$PC1 > 0))
outliers
DESeq2Table <- DESeq2Table[, !(colnames(DESeq2Table) %in% outliers)] 

#Reprise du heatmap et PCA sans les outliers
rld <- rlogTransformation(DESeq2Table, blind=TRUE)

distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)

rownames(mat) <-  colData(rld)$condition
colnames(mat) <-  colData(rld)$sampleNO

hmcol <- colorRampPalette(brewer.pal(9, "OrRd"))(255)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))

#PCA
Pvars <- rowVars(assay(rld))
select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, 
                                                      length(Pvars)))]

PCA <- prcomp(t(assay(rld)[select, ]), scale = F)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)

dataGG = data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], 
                    PC3 = PCA$x[,3], PC4 = PCA$x[,4], 
                    sampleNO = colData(rld)$sampleNO,
                    condition = colData(rld)$condition)

(qplot(PC1, PC2, data = dataGG, color =  condition, 
       main = "PC1 vs PC2, top variable genes, no outliers", size = I(6))
  + labs(x = paste0("PC1, VarExp:", round(percentVar[1],4)),
         y = paste0("PC2, VarExp:", round(percentVar[2],4)))
  + scale_colour_brewer(type="qual", palette=4)
)

##Lets stop here and try with my data 
sessionInfo()