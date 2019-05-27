##Test la visualisation des données de RNAseq
#Installation/loading des packages
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

##Load des données et mise en forme dans un objet DESeqDataSet
#12 fichiers sous la forme: IDEnsembl comptes, dans ~/R/PCA-cells/counts/
#Ouvrir les fichiers dans l'Environement de travail
library(tidyverse)
setwd('~/R/PCA-cells/counts')
temp <- list.files(pattern="*.counts")
subdata <- sub(pattern = "*.mdup.counts$", "", temp)
listdata <- list()
for (i in 1:length(temp))
{
  listdata[[i]]<-read_table2(temp[i], col_names = FALSE)
}
#Faire un full join sur la colonne X1 pour rassembler toutes les données dans la même table 
jdt <- reduce(listdata, full_join, by = "X1")
#retirer les 5 dernières lignes
jdt <- jdt[-c(63678:63682),]
#renommer les colonnes pour assurer le suivi des données
names(jdt)[2:13] = subdata
jdt.df <- as.data.frame(jdt)

##Retourne dans le wd d'origine 
setwd('~/R/PCA-cells/R-rnaseq')
#save la table jdt
write.table(jdt.df, '~/R/PCA-cells/R-rnaseq/20190527.deseqtest.cellLine.txt')

###Transformer la table en objet DESeq
head(jdt, 2)
#Transform la 1ere colonne en nom de lignes
jdt.df <- jdt.df %>% remove_rownames %>% column_to_rownames(var="X1")
jdt.mx <- as.matrix(jdt.df, row.names = 'X1') #transforme l'array en matrice
#Création de l'information de colonne
condition <- factor(c(rep("mut", 3), rep("WT", 3), rep("sh", 3), rep("ctl", 3)))
coldata <- data.frame(row.names=colnames(jdt.mx), condition)
#Vérifie que l'ordre de jdt.mx et coldata est le même
all(rownames(coldata) == colnames(jdt.mx)) #TRUE 
#Crée l'objet DESeqDateSet
dds <- DESeqDataSetFromMatrix(countData = jdt.mx,
                              colData = coldata,
                              design = ~ condition)
dds

#Pre-filtering du df, retire tous les gènes de moins de 10 reads
keep <- rowSums(counts(dds)) >= 10 #Totallement arbitraire
dds <- dds[keep, ]

#Établisssement d'un facteur de référence
dds$condition <- factor(dds$condition, 'ctl')

#Differential expression analysis 
dea.dds <- DESeq(dds)
resWT <- results(dea.dds, contrast = c('condition', 'WT', 'ctl'))
resmut <- results(dea.dds, contrast = c('condition', 'mut', 'ctl'))
resultsNames(dea.dds)

#Visualisation et ranking des resultats
#BiocManager::install('apeglm')
library(apeglm)
resLFC.WT <- lfcShrink(dea.dds, coef="condition_WT_vs_ctl", type="apeglm")
resLFC.mut <- lfcShrink(dea.dds, coef="condition_mut_vs_ctl", type="apeglm")
resLFC.WT
#Réarranger les résultats par pValue
resWTOrdered <- resWT[order(resWT$pvalue),]
summary(resWT)
#pValue plus petite de 0,1
sum(res$padj < 0.1, na.rm=TRUE)

#MA-plot for lfc results 
pdf("MA-plot_WT-mut.pdf")
par(mfrow=c(1,2), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC.WT, xlim=xlim, ylim=ylim, main="WildType") #rouge = pValue < 0,1
plotMA(resLFC.mut,xlim=xlim, ylim=ylim, main="Mutant")
dev.off()

#REpésentation d'un seul gène dans l'array, exemple celui avec la plus petite pValue
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

# Plot dispersions
png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dea.dds, main="Dispersion plot")
dev.off()


#Visualisation rapide des données (explorative)
vsd <- vst(dea.dds, blind=FALSE)
plotPCA(vsd)
rld <- rlog(dds, blind=FALSE)
plotPCA(rld)
#package supplémentaire
install.packages('pheatmap')
library(pheatmap)
#heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "OrRd")) )(255)
#Heatmap 
pdf("Heatmap_WT-mut.pdf")
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()
