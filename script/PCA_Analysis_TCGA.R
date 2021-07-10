
setwd('~/miRNomes/')

meta.tcga <- readRDS('shinyApp/data/Metadata_TCGA.RDS')
mir.tcga <- readRDS('shinyApp/data/miRNA_Expression_TCGA.RDS')

mir.annotation <- readRDS('shinyApp/data/miRBase_10.0_22.RDS')

###
projects <- names(meta.tcga)
projects

pca.tcga <- list()

for (prj in projects) {
  
  print (prj)
  
  expr <- mir.tcga[[prj]]
  meta <- meta.tcga[[prj]]

  keep <- rowSums(expr > 0) >= 0.5*ncol(expr)
  expr.med <- apply(expr[keep,], 1, median)
  
  o <- order(expr.med, decreasing = T)
  expr.med <- expr.med[o]
  
  mir.id <- names(expr.med)
  mir.name <- mir.annotation[mir.id,]$Name
  
  dataForPCA <- expr[mir.id,]
  
  filter <- which(rowSums(is.na(dataForPCA))>0)
  
  if (length(filter)>0) {
    dataForPCA <- t(scale(t(dataForPCA[-filter,])))
  } else {
    dataForPCA <- t(scale(t(dataForPCA)))
  }
  
  
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
                               Group=factor(meta$clinical_stage), # , levels=group.levels
                               Disease.Status=meta$sample_type,
                               stringsAsFactors = F)
  
  pca.tcga[[prj]] <- dataForPCAPlot
  
}


saveRDS(pca.tcga, file='PCA.Analysis.TCGA.RDS')
