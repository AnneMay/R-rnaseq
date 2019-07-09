##Filtres pour la "normalisation" des données du RNA-seq 
#Désire obtenir 3 listes: WT-Vecteur; WT-mut; WT-mut-corrigee
#Liste de départ: Résultats Deseq2 wt-vs-sramble et P904L-vs-wt

#  Chargement des library 
library(tidyverse)

#  Chargement des deux tables de types résultat DEseq2 
wt_vs_scramble_deseq2 <- read_table2("C:/Users/anmay/OneDrive - Universite de Montreal/Documents/Maîtrise/DNMT3A/Résultats/Résultats Lenti+sh/RNA-seq/WT_vs_Scramble/wt-vs-scramble.deseq2.tsv")
#P904L_vs_wt_deseq2 <- read_table2("C:/Users/anmay/OneDrive - Universite de Montreal/Documents/Maîtrise/DNMT3A/Résultats/Résultats Lenti+sh/RNA-seq/Fichier_RNA-seq/anne-marie_results/mut.deseq2.tsv")
mut_lenti_deseq2 <- read_delim("~/R/PCA-cells/mut-lenti.deseq2.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
wantedEnsemblID=P904L_vs_wt_deseq2[ ,1]
#Transform la 1ere colonne en nom de lignes
mut_lenti_deseq2.df <- mut_lenti_deseq2 %>% remove_rownames %>% column_to_rownames(var="gene_name")
mut_lenti_deseq2.dfc <- filter(mut_lenti_deseq2, gene_name %in% wantedEnsemblID$gene_name)

# Filtre 1: Tri des valeur de pValue ajusté < 0.05 (produit WT-vecteur et WT-mut)
#établissement du threshold 
wt_vs_scramble_deseq2$threshold = as.factor(wt_vs_scramble_deseq2$padj < 0.05)
#P904L_vs_wt_deseq2$threshold = as.factor(P904L_vs_wt_deseq2$padj < 0.05)
mut_lenti_deseq2.dfc$threshold = as.factor(mut_lenti_deseq2.dfc$padj < 0.05)

#retrait des valeurs FALSE
WT_vecteur <- filter(wt_vs_scramble_deseq2, threshold == TRUE)
#WT_mut <- filter(P904L_vs_wt_deseq2, threshold == TRUE)
Mut_vecteur <- filter(mut_lenti_deseq2.dfc, threshold == TRUE)

# Filtre 2: Impact du mutant seulement, output = 230 gènes
#outer join sur les symboles de gènes entre WT-Vecteur et WT-mut (WT-mut-corrigee partie A: corrA)
corrA <- as.data.frame(setdiff(WT_mut$gene_name, WT_vecteur$gene_name))
corrA <- corrA %>% rename("gene_name" = "setdiff(WT_mut$gene_name, WT_vecteur$gene_name)")
corrA <- corrA %>% left_join(WT_mut, by = "gene_name")

###  Gène différent entre Mutant_vecteur et WT_vecteur, pas vraiment ce qu'on veut, output = 873 gènes
# corrA2 <- as.data.frame(setdiff(Mut_vecteur$gene_name,WT_vecteur$gene_name))
# corrA2 <- corrA2 %>% rename("gene_name" = "setdiff(Mut_vecteur$gene_name, WT_vecteur$gene_name)")
# corrA2 <- corrA2 %>% left_join(WT_mut, by = "gene_name")

#test pour corrA2, output = 747 gènes
# test <- as.data.frame(setdiff(corrA2[,1], WT_mut$gene_name))
# test <- test %>% rename("gene_name" = "setdiff(corrA2[, 1], WT_mut$gene_name)")
# test <- test %>% left_join(WT_mut, by = "gene_name")

# Filtre 3: Gène subissant un impact mixte entre le mutant et le WT
#étape 1: établir la liste des gènes concernés: Inner join entre WT-Vecteur et WT-mut
intersect <- inner_join(WT_vecteur, WT_mut, by = "gene_name", suffix = c("1", "2"))
intersect2 <- inner_join(WT_vecteur, Mut_vecteur, by = "gene_name", suffix = c("1", "2"))
##Filtre 3 alternatif avec l'intersect 1, mais les valeurs de FC du Mut-plenti
wanted=intersect[ ,1]
Mut_vecteur.i <- filter(mut_lenti_deseq2.dfc, gene_name %in% wanted$gene_name)
intersect <- Mut_vecteur.i %>% left_join(WT_vecteur, by = "gene_name", suffix = c("1", "2")) #1= mut, 2= WT

#Renmaing and reorganising the element: gene_name, FC1 (vecteur), FC2 (mut)
intersect.r <- intersect %>% select(gene_name, "FC1" = log2FoldChange1, "FC2" = log2FoldChange2)
intersect.r2 <- intersect2 %>% select(gene_name, "FC1" = log2FoldChange1, "FC2" = log2FoldChange2)

#étape 2: Ressortir les statistiques de FC pour les deux comparaisons (WT-vecteur, WT-mut), on s'intéresse surtout au: 
#Range, quartile, moyenne.
resume <- summarise(intersect.r, avg = mean(cbind(intersect.r$FC1, intersect.r$FC2)), med_FC = median(cbind(intersect.r$FC1, intersect.r$FC2)),
                    min_FC = min(cbind(intersect.r$FC1, intersect.r$FC2)), max_FC = max(cbind(intersect.r$FC1, intersect.r$FC2)))
FC <- as.vector(as.matrix(intersect.r[,c("FC1", "FC2")]))
FCpos <- FC[FC > 0]
FCneg <- FC[FC < 0]

#Étape 3: Définir le threshold de FC 
FCthresholdpos <- quantile(FCpos, 0.05)
FCthresholdneg <- quantile(FCneg, 0.95)

#étape 4: Loop qui permet de tester chaque paire de FC dans le threshold
listA <- list()
listB <- list()
for(i in 1:nrow(intersect.r)){
  if ((as.numeric(intersect.r$FC1[i]) >= FCthresholdneg & as.numeric(intersect.r$FC2[i]) >= FCthresholdneg) || 
      (as.numeric(intersect.r$FC1[i]) <= FCthresholdpos & as.numeric(intersect.r$FC2[i]) <= FCthresholdpos)) {
    listA <- c(listA, intersect.r$gene_name[i])
  } else {listB <- c(listB, intersect.r$gene_name[i])}
}
# listC <- list()
# for(i in 1:nrow(intersect.r)){
#   if ((as.numeric(intersect.r$FC1[i]) >= 0 & as.numeric(intersect.r$FC2[i]) >= 0) || 
#       (as.numeric(intersect.r$FC1[i]) <= 0 & as.numeric(intersect.r$FC2[i]) <= 0)) {
#     listA <- c(listA, intersect.r$gene_name[i])
#   } else {listC <- c(listC, intersect.r$gene_name[i])}
# }
listD <- list()
for(i in 1:nrow(intersect.r2)){
  if ((as.numeric(intersect.r2$FC1[i]) >= FCthresholdneg & as.numeric(intersect.r2$FC2[i]) >= FCthresholdneg) ||
      (as.numeric(intersect.r2$FC1[i]) <= FCthresholdpos & as.numeric(intersect.r2$FC2[i]) <= FCthresholdpos)) {
    listA <- c(listA, intersect.r2$gene_name[i])
  } else {listD <- c(listD, intersect.r2$gene_name[i])}
}
#étape 5: réassocier les gènes avec leur valeur de FC (WT-mut) (si ça c'est perdu)
#(liste WT-mut-corrigee partie B: corrB)

listB <- as.data.frame(t(t(listB)))
listBc <- listB %>% rename("gene_name" = V1) %>% rownames_to_column() %>% select("gene_name") %>%  
  mutate("gene_name" = as.character("gene_name"))
corrB <- listBc %>% left_join(WT_mut, by = "gene_name")

### Combiner les deux listes de gènes: WT-mut-corrigee
WT_mut_corr <- bind_rows(corrA, corrB, ID = names("source"))

##Listes alternatives 
listD <- as.data.frame(t(t(listD)))
listDc <- listD %>% rename("gene_name" = V1) %>% rownames_to_column() %>% select("gene_name") %>%  
  mutate("gene_name" = as.character("gene_name"))
corrD <- listDc %>% left_join(WT_mut, by = "gene_name")
WT_mut_corr2 <- bind_rows(corrA2, corrD, ID = names("source"))

listD <- as.data.frame(t(t(listD)))
listDc <- listD %>% rename("gene_name" = V1) %>% rownames_to_column() %>% select("gene_name") %>%  
  mutate("gene_name" = as.character("gene_name"))
corrD <- listDc %>% left_join(WT_mut, by = "gene_name")
WT_mut_test <- bind_rows(test, corrD, ID = names("source"))
