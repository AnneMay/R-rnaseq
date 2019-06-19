#Test loop de CorrectionRNASeq.r 

#test dataset 
testDF <- matrix(c(sample(-100:100, 40, replace = FALSE)), ncol = 2, byrow = TRUE)
colnames(testDF) <- c( "FC1", "FC2")
rownames(testDF) <- c(LETTERS[seq( from = 1, to = 20 )])
testDF <- as.data.frame(testDF)
testDF <- testDF %>% rownames_to_column() 
testDF <- rename(testDF, "genesymbol" = rowname)


###Loop
#étape 2: Ressortir les statistiques de FC pour les deux comparaisons (WT-vecteur, WT-mut), on s'intéresse surtout au: 
#Range, quartile, moyenne.
FC <- as.vector(as.matrix(testDF[,c("FC1", "FC2")]))
FCpos <- FC[FC > 0]
FCneg <- FC[FC < 0]
#q2pos = quantile(FCpos, 0.05)
#q3neg = quantile(FCneg, 0.95)
resume <- summarise(testDF, avg = mean(cbind(testDF$FC1, testDF$FC2)), med_FC = median(cbind(testDF$FC1, testDF$FC2)),
                    min_FC = min(cbind(testDF$FC1, testDF$FC2)), max_FC = max(cbind(testDF$FC1, testDF$FC2)))


#Étape 3: Définir le threshold de FC 
FCthresholdpos <- quantile(FCpos, 0.05)
FCthresholdneg <- quantile(FCneg, 0.95)

#étape 4: Loop qui permet de tester chaque paire de FC dans le threshold
listA <- list()
listB <- list()
for(i in 1:nrow(testDF)){
  if ((as.numeric(testDF$FC1[i]) >= FCthresholdpos & as.numeric(testDF$FC2[i]) >= FCthresholdpos) || 
      (as.numeric(testDF$FC1[i]) <= FCthresholdneg & as.numeric(testDF$FC2[i]) <= FCthresholdneg)) {
    listA <- c(listA, testDF$genesymbol[i])
  } else {listB <- c(listB, testDF$genesymbol[i])}
}
#Étape 5: Retransformer les listes A et B en DF avec les informations de Deseq WT-Mut seulement
