
setwd('~/miRNomes/')

meta.tcga <- readRDS('shinyApp/data/Metadata_TCGA.RDS')
mir.tcga <- readRDS('shinyApp/data/miRNA_Expression_TCGA.RDS')

mir.annotation <- readRDS('shinyApp/data/miRBase_10.0_22.RDS')

###
projects <- names(meta.tcga)
projects

roc.list <- list()

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
  
  dataForROCAnalysis <- expr[mir.id,]
  
  filter <- which(rowSums(is.na(dataForROCAnalysis))>0)
  
  if (length(filter)>0) {
    dataForROCAnalysis <- dataForROCAnalysis[-filter,]
  }
  
  group <- meta$sample_type
  
  dataForForestPlot <- c()
  
  if (length(unique(group))==2) {
    
    for (mir in mir.id) {
      
      roc.test <- roc(group, dataForROCAnalysis[mir,], plot=FALSE, ci=TRUE, auc=TRUE)
      ci.auc <- roc.test$ci
      
      auc <- ci.auc[2]
      auc.ci.lower95 <- ci.auc[1]
      auc.ci.upper95 <- ci.auc[3]
      
      auc <- format(auc, digits = 2, nsmall=2)
      auc.ci.lower95 <- format(auc.ci.lower95, digits = 2, nsmall=2)
      auc.ci.upper95 <- format(auc.ci.upper95, digits = 2, nsmall=2)
      
      
      pred <- prediction(dataForROCAnalysis[mir,], group)
      perf <- performance(pred, measure = "tpr", x.measure = "fpr")
      
      FPR <- perf@x.values[[1]]
      TPR <- perf@y.values[[1]]
      
      #df <- data.frame(FPR,TPR)
      
      auc.test <- wilcox.test(FPR, TPR, alternative = 'two.sided')
      pvalue <- formatC(auc.test$p.value, format = 'e', digits = 2)
      
      dataForForestPlot <- rbind(dataForForestPlot, 
                                 c(auc, auc.ci.lower95, auc.ci.upper95, pvalue))
      
    }
    
    dataForForestPlot <- apply(dataForForestPlot, 2, as.numeric)
    
    dataForForestPlot <- data.frame(miRNA.Accession=mir.id,
                                    miRNA.ID=mir.name,
                                    dataForForestPlot,
                                    row.names = NULL,
                                    stringsAsFactors = F)
    
    colnames(dataForForestPlot)[3:6] <- c('AUC','Lower95','Upper95','P.Value')
    
    o <- order(dataForForestPlot$P.Value, decreasing = F)
    dataForForestPlot <- dataForForestPlot[o,]
    
    o <- order(dataForForestPlot$AUC, decreasing = T)
    dataForForestPlot <- dataForForestPlot[o,]

    
  } else {
    dataForForestPlot <- data.frame(matrix(ncol = 5, nrow = 0))
    dataForForestPlot
    colnames(dataForForestPlot) <- c('AUC','Lower95','Upper95','P.Value','miRNA.ID')
    
  }
  
  roc.list[[prj]] <- dataForForestPlot
}

saveRDS(roc.list, file='ROC.Analysis.TCGA.RDS')
