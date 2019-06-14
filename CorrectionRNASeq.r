##Filtres pour la "normalisation" des données du RNA-seq 
#Désire obtenir 3 listes: WT-Vecteur; WT-mut; WT-mut-corrigee
#Liste de départ: Résultats Deseq2 wt-vs-sramble et P904L-vs-wt

#  Chargement des library 
library(tidyverse)

#  Chargement des deux tables de types résultat DEseq2 
wt_vs_scramble_deseq2 <- read_table2("C:/Users/anmay/OneDrive - Universite de Montreal/Documents/Maîtrise/DNMT3A/Résultats/Résultats Lenti+sh/RNA-seq/WT_vs_Scramble/wt-vs-scramble.deseq2.tsv")
P904L_vs_wt_deseq2 <- read_table2("C:/Users/anmay/OneDrive - Universite de Montreal/Documents/Maîtrise/DNMT3A/Résultats/Résultats Lenti+sh/RNA-seq/P904L_vs_wt/P904L_vs_wt.deseq2.tsv")

# Filtre 1: Tri des valeur de pValue ajusté < 0.05 (produit WT-vecteur et WT-mut)
#établissement du threshold 
wt_vs_scramble_deseq2$threshold = as.factor(wt_vs_scramble_deseq2$padjst < 0,05)
P904L_vs_wt_deseq2$threshold = as.factor(P904L_vs_wt_deseq2$padjst < 0,05)

#retrait des valeurs FALSE
WT_vecteur <- filter(wt_vs_scramble_deseq2, threshold == TRUE)
WT_mut <- filter(P904L_vs_wt_deseq2, threshold == TRUE)

# Filtre 2: Impact du mutant seulement
#outer join sur les symboles de gènes entre WT-Vecteur et WT-mut (WT-mut-corrigee partie A: corrA)
corrA <- setdiff(WT_mut, WT_vecteur)

# Filtre 3: Gène subissant un impact mixte entre le mutant et le WT
#étape 1: établir la liste des gènes concernés: Inner join entre WT-Vecteur et WT-mut
intercept <- inner_join(WT_vecteur, WT_mut, by = "gene symbol", suffix = c("1", "2"))

#étape 2: Ressortir les statistiques de FC pour les deux comparaisons (WT-vecteur, WT-mut), on s'intéresse surtout au: 
#Range, quartile, moyenne.
summarise(intercept)

#Étape 3: Définir le threshold de FC 
FCthresholdpos <- as_factor()
FCthresholdneg <- as_factor()

#étape 4: Loop qui permet de tester chaque paire de FC dans le threshold
# listA <- list()
# listB <- list()
# for i in range(intercept[ , ]){
#   if (intercept$FC1[i] >= FCthresholdpos && intercept$FC2[i] >= FCthresholdpos) || 
#   (intercept$FC1[i] <= FCthresholdneg && intercept$FC2[i] <= FCthresholdneg) {
#     print i in listA
#   } else {print i in listB}
# }
 
#étape 5: réassocier les gènes avec leur valeur de FC (WT-mut) (si ça c'est perdu)
#(liste WT-mut-corrigee partie B: corrB)
listB %>% left_join(intercept, by = "gene symbol")
corrB <- select(listB, #nom colone sauf FC WT-vecteur
                )

### Combiner les deux listes de gènes: WT-mut-corrigee
WT_mut_corr <- bind_rows(corrA, corrB, ID = "source")