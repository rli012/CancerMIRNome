
############################################

setwd('~/Publications/CancerMIRNome/')

meta.tcga <- readRDS('shinyApp/data/Metadata_TCGA.RDS')
mir.tcga <- readRDS('shinyApp/data/miRNA_Expression_TCGA.RDS')

###
projects <- names(meta.tcga)
projects

View(meta.tcga[['TCGA-COAD']])

rownames(meta.tcga[['TCGA-COAD']]) == colnames(mir.tcga[['TCGA-COAD']])
mir.tcga[['TCGA-COAD']][1:5,1:5]


for (prj in projects) {
  
  print (prj)
  
  expr <- mir.tcga[[prj]]
  meta <- meta.tcga[[prj]]
  
  eSet <- ExpressionSet(assayData = as.matrix(expr),
                        phenoData = AnnotatedDataFrame(meta))
  
  saveRDS(eSet, file=paste0('~/Publications/OncomiRNomeDB/www/downloads/TCGA/', prj, '_eSet.RDS'))

}


prad <- readRDS('~/Publications/OncomiRNomeDB/www/downloads/TCGA-PRAD_eSet.RDS')
exprs(prad)[1:5,1:5]
pData(prad)[1:5,1:5]

dim(exprs(prad))

all(rownames(pData(prad))==colnames(exprs(prad)))


##########

ccma.primary <- readRDS('shinyApp/data/miRNomes_Datasets_Primary.RDS')
ccma.primary

nrow(ccma.primary)

expr.ccma <- readRDS('shinyApp/data/miRNomes_Expression.RDS')
meta.ccma <- readRDS('shinyApp/data/miRNomes_Metadata.RDS')

names(expr.ccma) == names(meta.ccma)
names(expr.ccma)

meta.ccma[['GSE122497']]

for (idx in 1:nrow(ccma.primary)) {
  dataset <- as.character(ccma.primary[idx,'Dataset'])
  
  print (dataset)
  
  meta <- meta.ccma[[dataset]]
  expr <- expr.ccma[[dataset]]
  
  eSet <- ExpressionSet(assayData = as.matrix(expr),
                        phenoData = AnnotatedDataFrame(meta))
  
  saveRDS(eSet, file=paste0('~/Publications/OncomiRNomeDB/www/downloads/Circulating/', dataset, '_eSet.RDS'))
  
  
}

dim(expr)

###
mir.annotation <- readRDS('shinyApp/data/miRBase_10.0_22.RDS')
View(mir.annotation)

saveRDS(mir.annotation, '~/Publications/OncomiRNomeDB/www/downloads/Annotation/miRNA_Annotation_miRBase_Release10.0_to_Release22.RDS')

