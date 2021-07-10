
###
setwd('~/CancerMIRNome/')

meta.tcga <- readRDS('shinyApp/data/Metadata_TCGA.RDS')
mir.tcga <- readRDS('shinyApp/data/miRNA_Expression_TCGA.RDS')

mir.annotation <- readRDS('shinyApp/data/miRBase_10.0_22.RDS')

###
projects <- names(meta.tcga)
projects

feature.plot.list <- list()
feature.table.list <- list()


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
  
  if (length(unique(group))==2 & (! prj %in% c('TCGA-SKCM','TCGA-THYM','TCGA-CESC'))) {
    print (! prj %in% c('TCGA-SKCM','TCGA-THYM','TCGA-CESC'))
    
    cvfit<-cv.glmnet(x=t(dataForROCAnalysis),y=group, alpha=1, 
                     family = "binomial", type.measure="class") #class
    
    feature.plot.list[[prj]] <- cvfit
    
    coef.min<-coef(cvfit,s="lambda.min")
    feature <- data.frame(coef.min@Dimnames[[1]],matrix(coef.min), stringsAsFactors = F)
    colnames(feature) <- c('miRNA.Accession','Coefficient')
    
    keep <- which(feature$Coefficient!=0)
    feature <- feature[keep,]
    
    feature <- data.frame(miRNA.Accession=feature$miRNA.Accession,
                          miRNA.ID=mir.annotation[feature$miRNA.Accession,]$Name,
                          Coefficient=feature$Coefficient,
                          row.names = NULL,
                          stringsAsFactors = F)
    
    o <- order(abs(feature$Coefficient), decreasing = T)
    feature <- feature[o,]
    

  } else {
    
    anno <- data.frame(x = 1, y = 1, label = 'Warning: no feature is selected')
    
    p <- ggplot() + geom_text(data=anno, aes(x,y, label=label), 
                              color='red', size=5.5, fontface='bold.italic') +
      xlim(0,2) + ylim(0,2) +
      theme_bw()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank())
    
    feature.plot.list[[prj]] <- p
    
    feature <- data.frame(matrix(ncol = 3, nrow = 0))
    colnames(feature) <- c('miRNA.Accession','miRNA.ID','Coefficient')
    
    
  }
  
  feature.table.list[[prj]] <- feature
  
}

saveRDS(feature.plot.list, file='Lasso.Feature.Plot.RDS')
saveRDS(feature.table.list, file='Lasso.Feature.Table.RDS')



plot(feature.plot.list[[1]])

feature.plot.list[[2]]



