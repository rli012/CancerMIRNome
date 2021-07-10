
setwd('~/miRNomes/')

meta.tcga <- readRDS('shinyApp/data/Metadata_TCGA.RDS')
mir.tcga <- readRDS('shinyApp/data/miRNA_Expression_TCGA.RDS')

mir.annotation <- readRDS('shinyApp/data/miRBase_10.0_22.RDS')

###
projects <- names(meta.tcga)
projects

expr.high.list <- list()

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
  
  expr.high <- data.frame(miRNA.Accession=mir.id, miRNA.ID=mir.name, 
                          Median=expr.med, 
                          Rank=1:length(mir.id),
                          stringsAsFactors = F)
  
  expr.high.list[[prj]] <- expr.high
  
}

expr.high.list[[prj]]  

saveRDS(expr.high.list, file='Highly.Expressed.miRNAs.TCGA.RDS')
