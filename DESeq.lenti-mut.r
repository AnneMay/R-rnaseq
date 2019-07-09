###Test avec la comparaison mut-vecteur
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

#définir le fichier design et data
data <- select(jdt.df, X1, ACHN_DNMT3AmutP904L_1, ACHN_DNMT3AmutP904L_2, ACHN_DNMT3AmutP904L_3, 
               ACHN_Shscramble_1, ACHN_Shscramble_2, ACHN_Shscramble_3)
design <- read_delim("~/R/PCA-cells/mut-vs-scramble.design.txt", 
                     "\t", escape_double = FALSE, trim_ws = TRUE)

#Manque l'étape de retrait des ID gène non désiré, il me manque un fichier 


design$group=as.factor(design$group)
#Transform la 1ere colonne en nom de lignes
jdt.df <- data %>% remove_rownames %>% column_to_rownames(var="X1")
data <- as.matrix(jdt.df, row.names = 'X1') #transforme l'array en matrice


dds <- DESeqDataSetFromMatrix(countData = data, colData = design, design = ~ group)
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),]

d1=cbind(rownames(resOrdered),data.frame(resOrdered))
colnames(d1)[1]=c("gene_name")
head(counts(dds, normalized=T))
d2=cbind(rownames(counts(dds, normalized=T)),data.frame(counts(dds, normalized=T)))
colnames(d2)[1]=c("gene_name")
colnames(d2)[2:ncol(d2)]=colnames(data)
d3=merge(d1,d2,by.x="gene_name",by.y="gene_name")
d3[8:ncol(d3)]=round(d3[8:ncol(d3)],2)
d3=d3[ order(d3$padj), ]

# Merge with gene function/GO term

queryVec=d3$gene_name
names <- queryMany(queryVec, scopes=c("ensembl.gene"), fields=c("symbol", "name"), species='human')
sym_names <- names[,c('query','symbol','name')]
colnames(sym_names) <- c('gene_name','symbol','description')
go <- queryMany(queryVec, scopes='ensembl.gene', fields="go.MF", species='human')
go_data = NULL
for (i in c(1:nrow(go))){
  g <- go[i,'query']
  if(length(go[i,'go.MF'][[1]]$term) >= 1){
    go_term = go[i,'go.MF'][[1]]$term[1]
    go_data <- rbind(go_data,data.frame(gene_name=g,term=go_term))
  }
}

final=merge(sym_names,go_data,all=T)
final=merge(d3,final,by.x=1,by.y=1)
final=subset(final, !duplicated(final$gene_name))
# final=d3
final=final[ order(final$padj), ]
final$baseMean=round(final$baseMean,3)
final$log2FoldChange=round(final$log2FoldChange,3)
final$lfcSE=round(final$lfcSE,3)
final$stat=round(final$stat,3)

write.table(final,"mut-lenti.deseq2.tsv", quote=F, row.names=F,sep="\t")
rld <- rlog(dds)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat)
hc <- hclust(distsRL)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)


pdf("RNASeq_Clustering_mut-lenti_June_2017.pdf")
heatmap.2(mat, Rowv=as.dendrogram(hc),symm=TRUE, trace="none",col = rev(hmcol), margin=c(13, 13), cexRow = 0.7, cexCol = 0.7)
plotPCA(rld, intgroup=c("group"))
dev.off()
