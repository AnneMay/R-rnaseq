# devtools::install_github("GuangchuangYu/clusterProfiler")
# devtools::install_github("GuangchuangYu/DOSE")
# devtools::install_github("GuangchuangYu/enrichplot")


library("org.Hs.eg.db")
library("ggplot2")
library(clusterProfiler)
library("DOSE")
library("pathview")


#Read deseq2 data
# deseqRes <- read.table("N91.deseq2.tsv",sep="\t",header=T,quote="")
deseqRes <- as.data.frame(WT_vecteur)
deseqRes$gene_name <- as.character(deseqRes$gene_name)
deseqRes <- deseqRes[grep("ENSG",deseqRes$gene_name),]


#Filter for padjusted value and extract appropriate column for ranking
deseqResFil <- deseqRes[!is.na(deseqRes$padj),]

#log2FoldChange
DEgenes <- deseqResFil[deseqResFil$padj<0.05,3]
names(DEgenes) <- deseqResFil[deseqResFil$padj<0.05,1]
DEgenes <- DEgenes[order(-DEgenes)]

#stat
DEgenesStat <- deseqResFil[deseqResFil$padj<0.05 & abs(deseqResFil$log2FoldChange)>1 ,"stat"]
names(DEgenesStat) <- deseqResFil[deseqResFil$padj<0.05 & abs(deseqResFil$log2FoldChange)>1,1]
DEgenesStat <- DEgenesStat[order(-DEgenesStat)]


#For gene enrichment
genes=as.integer(deseqResFil$padj<0.05)

names(genes)=deseqResFil$gene_name
genes[is.na(genes)]<-0

#Fetch Entrez ID
geneNamesOrg <- mapIds(org.Hs.eg.db,
                    keys=names(DEgenes), 
                    column="ENTREZID",
                    keytype="ENSEMBL",
                    multiVals="first")

geneNamesOrgStat <- mapIds(org.Hs.eg.db,
                       keys=names(DEgenesStat), 
                       column="ENTREZID",
                       keytype="ENSEMBL",
                       multiVals="first")
#



#Assign entrezID
DEgenesEntrez <- DEgenes
DEgenesEntrezStat <- DEgenesStat

names(DEgenesEntrez) <- geneNamesOrg
names(DEgenesEntrezStat) <- geneNamesOrgStat

DEgenesEntrez <- DEgenesEntrez[!is.na(names(DEgenesEntrez))]
DEgenesEntrez <- DEgenesEntrez[order(-DEgenesEntrez)]

DEgenesEntrezStat <- DEgenesEntrezStat[!is.na(names(DEgenesEntrezStat))]

#Fetch symbols
geneSymbolOrg <- mapIds(org.Hs.eg.db,
                       keys=names(DEgenes), 
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first")

DEgenesSymbolLfc <- DEgenes
names(DEgenesSymbolLfc) <- geneSymbolOrg


geneSymbolOrgStat <- mapIds(org.Hs.eg.db,
                        keys=names(DEgenesStat), 
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")

DEgenesSymbolStat <- DEgenesStat
names(DEgenesSymbolStat) <- geneSymbolOrgStat


#Run erichment analysis for GO stat
ego2Stat <- enrichGO(gene         = names(DEgenesStat),
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)


#Run gene set enrichment analysis for GO stat
ego4Stat <- gseGO(geneList     = DEgenesEntrezStat,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              nPerm        = 1000,
              minGSSize    = 15,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE,
              seed         = 223456)

####GMTs for other pathways

#Read GMT
#  hallm <- read.gmt("h.all.v6.2.entrez.gmt")
#  c6 <- read.gmt("c6.all.v6.2.entrez.gmt")

#Run GSEA on GMTs
#hallGseaStat <- GSEA(DEgenesEntrezStat, 
#                 TERM2GENE=hallm, 
#                 minGSSize = 10,
#                 maxGSSize = 500,
#                 pvalueCutoff = 0.05,
#                 verbose=F,
#                 seed = 223456)


#c6GseaStat <- GSEA(DEgenesEntrezStat, 
#                 TERM2GENE=c6, 
#                 minGSSize = 10,
#                 maxGSSize = 500,
#                 pvalueCutoff = 0.05,
#                 verbose=F,
#                 seed = 223456)





#Write GSEA results to file
egoxStat <- setReadable(ego4Stat, 'org.Hs.eg.db', 'auto')
#hallGseaStat <- setReadable(hallGseaStat, 'org.Hs.eg.db', 'ENTREZID')
#c6GseaStat <- setReadable(c6GseaStat, 'org.Hs.eg.db', 'ENTREZID')
write.csv(egoxStat,file="GO_GSEA_Stat.csv")
#write.csv(hallGseaStat,file="hall_GSEA_Stat.csv")
#write.csv(c6GseaStat,file="c6_GSEA_Stat.csv")


#Write GSEA plots to pdf
pdf(file="GSEA_GO_BP_Stat_analysis.pdf",onefile = TRUE,width = 18,height=8)
dp5 <- dotplot(ego4Stat, showCategory=30) + ggtitle("dotplot for GO")
#p2<-barplot(ego4P, showCategory=30)
p3<-cnetplot(egoxStat, foldChange=DEgenesSymbolLfc)
p4<-heatplot(egoxStat, foldChange=DEgenesSymbolLfc)
p5<-emapplot(ego4Stat)
print(dp5)
#print(p2)
print(p3)
print(p4)
print(p5)
dev.off()


#hallGseaStat <- simplify(hallGseaStat)

#pdf(file="GSEA_hallmark_analysis_Stat.pdf",onefile = TRUE,width = 18,height=8)
#dp5 <- dotplot(hallGseaStat, showCategory=30) + ggtitle("dotplot for hallmark")
#p3<-cnetplot(hallGseaStat, foldChange=DEgenesSymbolLfc)
#p4<-heatplot(hallGseaStat, foldChange=DEgenesSymbolLfc)
#p5<-emapplot(hallGseaStat)
#print(dp5)
#print(p3)
#print(p4)
#print(p5)
#dev.off()


#c6GseaStat <- simplify(c6GseaStat)
#pdf(file="GSEA_c6_analysis_Stat.pdf",onefile = TRUE,width = 18,height=8)
#dp5 <- dotplot(c6GseaStat, showCategory=30) + ggtitle("dotplot for c6")
#p3<-cnetplot(c6GseaStat, foldChange=DEgenesSymbolLfc)
#p4<-heatplot(c6GseaStat, foldChange=DEgenesSymbolLfc)
#p5<-emapplot(c6GseaStat)
#print(dp5)
#print(p3)
#print(p4)
#print(p5)
# dev.off()




#kk <- enrichKEGG(gene = names(DEgenesEntrez), organism = 'hsa', pvalueCutoff = 0.05)
kk <- gseKEGG(geneList     = DEgenesEntrezStat,
        organism     = 'hsa',
        nPerm        = 1000,
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 0.05,
        verbose      = FALSE)

kpath <- head(kk)[[1]][[1]]
write.csv(as.data.frame(kk),file=paste(currexp,".KEGG_OverRepresentationAnalysis.csv",sep=""))
for(i in as.data.frame(kk)[[1]]){
  
  pv <- pathview(gene.data  = DEgenesEntrez, pathway.id = i, species    = "hsa", limit      = list(gene=max(abs(DEgenesEntrez)), cpd=1),out.suffix="KEGG_pathview")
  file.remove(c(paste(i,".png",sep=""),paste(i,".xml",sep="")))
  
}


