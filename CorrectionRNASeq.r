##Filtres pour la "normalisation" des données du RNA-seq 
#Désire obtenir 3 listes: WT-Vecteur; WT-mut; WT-mut-corrigee
#Liste de départ: Résultats Deseq2 wt-vs-sramble et P904L-vs-wt

#  Chargement des library 
library(tidyverse)

#  Chargement des deux tables de types résultat DEseq2 
wt_vs_scramble_deseq2 <- read_table2("C:/Users/anmay/OneDrive - Universite de Montreal/Documents/Maîtrise/DNMT3A/Résultats/Résultats Lenti+sh/RNA-seq/WT_vs_Scramble/wt-vs-scramble.deseq2.tsv")
P904L_vs_wt_deseq2 <- read_table2("C:/Users/anmay/OneDrive - Universite de Montreal/Documents/Maîtrise/DNMT3A/Résultats/Résultats Lenti+sh/RNA-seq/Fichier_RNA-seq/anne-marie_results/mut.deseq2.tsv")

# Filtre 1: Tri des valeur de pValue ajusté < 0.05 (produit WT-vecteur et WT-mut)
#établissement du threshold 
wt_vs_scramble_deseq2$threshold = as.factor(wt_vs_scramble_deseq2$padj < 0.05)
P904L_vs_wt_deseq2$threshold = as.factor(P904L_vs_wt_deseq2$padj < 0.05)

#retrait des valeurs FALSE
WT_vecteur <- filter(wt_vs_scramble_deseq2, threshold == TRUE)
WT_mut <- filter(P904L_vs_wt_deseq2, threshold == TRUE)

# Filtre 2: Impact du mutant seulement
#outer join sur les symboles de gènes entre WT-Vecteur et WT-mut (WT-mut-corrigee partie A: corrA)
corrA <- as.data.frame(setdiff(WT_mut$gene_name, WT_vecteur$gene_name))
corrA <- corrA %>% rename("gene_name" = "setdiff(WT_mut$gene_name, WT_vecteur$gene_name)")
corrA <- corrA %>% left_join(WT_mut, by = "gene_name")

# Filtre 3: Gène subissant un impact mixte entre le mutant et le WT
#étape 1: établir la liste des gènes concernés: Inner join entre WT-Vecteur et WT-mut
intersect <- inner_join(WT_vecteur, WT_mut, by = "gene_name", suffix = c("1", "2"))
#Renmaing and reorganising the element: gene_name, FC1 (vecteur), FC2 (mut)
intersect.r <- intersect %>% select(gene_name, "FC1" = log2FoldChange1, "FC2" = log2FoldChange2)

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
  if ((as.numeric(intersect.r$FC1[i]) >= FCthresholdpos & as.numeric(intersect.r$FC2[i]) >= FCthresholdpos) || 
      (as.numeric(intersect.r$FC1[i]) <= FCthresholdneg & as.numeric(intersect.r$FC2[i]) <= FCthresholdneg)) {
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
#étape 5: réassocier les gènes avec leur valeur de FC (WT-mut) (si ça c'est perdu)
#(liste WT-mut-corrigee partie B: corrB)

listB <- as.data.frame(t(listB)) 
listB <- listB %>% rename("gene_name" = V1) %>% rownames_to_column() %>% select(gene_name)
corrB <- listB %>% left_join(WT_mut, by = "gene_name")

### Combiner les deux listes de gènes: WT-mut-corrigee
WT_mut_corr <- bind_rows(corrA, corrB, ID = names("source"))