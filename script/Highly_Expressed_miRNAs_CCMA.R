
setwd('~/miRNomes/')

mir.annotation <- readRDS('shinyApp/data/miRBase_10.0_22.RDS')

ccma.primary <- readRDS('shinyApp/data/miRNomes_Datasets_Primary.RDS')
ccma.primary

nrow(ccma.primary)

expr.ccma <- readRDS('shinyApp/data/miRNomes_Expression.RDS')
meta.ccma <- readRDS('shinyApp/data/miRNomes_Metadata.RDS')


high.expr.ccma <- list()

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
  
  expr.high <- data.frame(miRNA.Accession=mir.id, miRNA.ID=mir.name, 
                          Median=expr.med, 
                          Rank=1:length(mir.id),
                          stringsAsFactors = F)
  
  high.expr.ccma[[dataset]] <- expr.high
  
  
}

saveRDS(high.expr.ccma, file='Highly.Expressed.miRNAs.CCMA.RDS')



