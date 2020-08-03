
meta.tcga <- readRDS('shinyApp/data/Metadata_TCGA.RDS')
mir.tcga <- readRDS('shinyApp/data/miRNA_Expression_TCGA.RDS')


projects <- names(meta.tcga)
projects


coxph.list <- list()

for (prj in projects) {
  
  print (prj)

  expr <- mir.tcga[[prj]]
  meta <- meta.tcga[[prj]]
  
  
  keep <- rowSums(expr > 0) >= 0.5*ncol(expr)
  
  samples <- which(meta$sample_type=='Tumor')
  
  expr <- expr[keep,samples]
  meta <- meta[samples,]

  mir.id <- rownames(expr)
  mir.name <- mir.annotation[mir.id,]$Name

  
  
  
  
  
  
  
    
}


