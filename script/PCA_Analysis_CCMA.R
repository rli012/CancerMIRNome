
setwd('~/miRNomes/')

mir.annotation <- readRDS('shinyApp/data/miRBase_10.0_22.RDS')

ccma.primary <- readRDS('shinyApp/data/miRNomes_Datasets_Primary.RDS')
ccma.primary

expr.ccma <- readRDS('shinyApp/data/miRNomes_Expression.RDS')
meta.ccma <- readRDS('shinyApp/data/miRNomes_Metadata.RDS')


pca.ccma <- list()

for (idx in 1:nrow(ccma.primary)) {
  dataset <- as.character(ccma.primary[idx,'Dataset'])
  
  print (dataset)
  
  meta <- meta.ccma[[dataset]]
  expr <- expr.ccma[[dataset]]

  expr.med <- apply(expr, 1, median)
  
  filter <- which(!names(expr.med) %in% rownames(mir.annotation))
  
  if (length(filter)>0) {
    expr.med <- expr.med[-filter]
  }
  
  o <- order(expr.med, decreasing = T)
  expr.med <- expr.med[o][1:500]
  
  mir.id <- names(expr.med)
  mir.name <- mir.annotation[mir.id,]$Name
  
  dataForPCA <- expr[mir.id,]
  
  filter <- which(rowSums(is.na(dataForPCA))>0)
  
  if (length(filter)>0) {
    dataForPCA <- dataForPCA[-filter,]
  }
  
  expr.sd <- apply(dataForPCA, 1, sd)
  filter <- which(expr.sd==0)
  
  if (length(filter)>0) {
    dataForPCA <- dataForPCA[-filter,]
  }
  
  dataForPCA <- t(scale(t(dataForPCA)))
  
  pcaResults <- prcomp(dataForPCA)
  sumpca <- summary(pcaResults)
  
  pc1 <- round(sumpca$importance[2,1]*100,2)
  pc2 <- round(sumpca$importance[2,2]*100,2)
  pc3 <- round(sumpca$importance[2,3]*100,2)
  
  dataForPCAPlot <- data.frame(PC1=pcaResults$rotation[,1],
                               PC2=pcaResults$rotation[,2],
                               PC3=pcaResults$rotation[,3],
                               pc1, pc2, pc3,
                               Sample=rownames(pcaResults$rotation),
                               Group=factor(meta$Group), # , levels=group.levels
                               Disease.Status=meta$Disease.Status,
                               stringsAsFactors = F)
  
  pca.ccma[[dataset]] <- dataForPCAPlot
  
}

saveRDS(pca.ccma, file='PCA.Analysis.CCMA.RDS')

