
################################## Data #####################################

# library(shiny)
# library(shinydashboard)
# library(shinydashboardPlus)
# library(shinyjs)
# #library(shinyWidgets)
# library(shinycssloaders)
# #library(readxl)
# library(ggplot2)
# library(Matrix)
# library(stringr)
# library(Biobase)
# library(survival)
# library(survminer)
# library(limma)
# library(edgeR)
# library(DT)
# library(dplyr)
# library(pROC)
# library(ROCR)
# library(digest)
# library(pheatmap)
# library(plotly)
# library(glmnet)
# library(heatmaply)
# library(htmlwidgets)
# library(slickR)
# 
# library(dashboardthemes)
# library(shinythemes)
# library(shinyalert)


### TCGA Datasets
tcga.datasets <- readRDS('data/TCGA_Projects.RDS')
idx <- which(tcga.datasets$Project=='TCGA-GBM')
tcga.datasets <- tcga.datasets[-idx,]


projects.tcga <- c("TCGA-ACC","TCGA-BLCA","TCGA-BRCA","TCGA-CESC",
                   "TCGA-CHOL","TCGA-COAD","TCGA-DLBC","TCGA-ESCA",#"TCGA-GBM",
                   "TCGA-HNSC","TCGA-KICH","TCGA-KIRC","TCGA-KIRP",
                   "TCGA-LAML","TCGA-LGG","TCGA-LIHC","TCGA-LUAD",
                   "TCGA-LUSC","TCGA-MESO","TCGA-OV","TCGA-PAAD",
                   "TCGA-PCPG","TCGA-PRAD","TCGA-READ","TCGA-SARC",
                   "TCGA-SKCM","TCGA-STAD","TCGA-TGCT","TCGA-THCA",
                   "TCGA-THYM","TCGA-UCEC","TCGA-UCS","TCGA-UVM")

projects.tcga.sub <- c("TCGA-BLCA","TCGA-BRCA","TCGA-CESC","TCGA-CHOL",
                       "TCGA-COAD","TCGA-ESCA","TCGA-HNSC","TCGA-KICH",
                       "TCGA-KIRC","TCGA-KIRP","TCGA-LIHC","TCGA-LUAD",
                       "TCGA-LUSC","TCGA-PAAD","TCGA-PCPG","TCGA-PRAD",
                       "TCGA-READ","TCGA-SKCM","TCGA-STAD","TCGA-THCA",
                       "TCGA-THYM","TCGA-UCEC")

### CCMA Datasets
ccma.datasets <- readRDS('data/miRNomes_Datasets.RDS')
ccma.datasets$Name <- paste0(ccma.datasets$Dataset, ': ', ccma.datasets$Title)
ccma.primary <- readRDS('data/miRNomes_Datasets_Primary.RDS')
ccma.primary <- ccma.primary[,c(1:3,5)]

### TCGA Data
meta.tcga <- readRDS('data/Metadata_TCGA.RDS')
mir.tcga <- readRDS('data/miRNA_Expression_TCGA.RDS')
#rna.tcga <- readRDS('data/RNAseq_Expression_TCGA.miRTarBase.RDS')

### TCGA Data Analysis
expr.high.tcga <- readRDS(file='data/Highly.Expressed.miRNAs.TCGA.RDS')

km.tcga <- readRDS('data/Survival.KM.TCGA.RDS')
coxph.tcga <- readRDS('data/Survival.CoxPH.TCGA.RDS')

lasso.tcga <- readRDS('data/Survival.Lasso.Feature.TCGA.RDS')
lasso.plot.tcga <- readRDS('data/Survival.Lasso.Plot.TCGA.RDS')

#risk.km.plot.tcga <- readRDS(file='data/Survival.KM.Risk.Plot.TCGA.RDS')
risk.km.data.tcga <- readRDS(file='data/Survival.KM.Risk.Data.TCGA.RDS')

tcga.feature.table <- readRDS('data/Lasso.Feature.Table.RDS')
tcga.feature.plot <- readRDS('data/Lasso.Feature.Plot.RDS')

roc.tcga <- readRDS('data/ROC.Analysis.TCGA.RDS')
surv.roc.plot.tcga <- readRDS(file='data/Survival.ROC.Risk.Plot.TCGA.RDS')

pca.tcga <- readRDS(file='data/PCA.Analysis.TCGA.RDS')

### Correlation/Functional Analysis
#cor.table <- readRDS('data/Correlation.miRTarBase.RDS')
enrichment.table <- readRDS('data/Enrichment.miRTarBase.RDS')



### CCMA Data Analysis
expr.high.ccma <- readRDS(file='data/Highly.Expressed.miRNAs.CCMA.RDS')

#expr.ccma <- readRDS('data/miRNomes_Expression.RDS')
meta.ccma <- readRDS('data/miRNomes_Metadata.RDS')

pca.ccma <- readRDS(file='data/PCA.Analysis.CCMA.RDS')




###
# meta.tcga <- readRDS('~/Documents/Publications/CancerMIRNomeAWS/data/Metadata_TCGA.RDS')
# meta <- do.call(rbind, meta.tcga)
# 
# sum(unlist(lapply(meta.tcga, nrow)))
# 
# length(meta.tcga)
# o <- order(names(meta.tcga))
# unlist(lapply(meta.tcga, nrow))[o]
# 
# saveRDS(meta, file='~/Documents/Publications/CancerMIRNomeAWS/data/Integrated_Sample_Metadata_TCGA.RDS')

