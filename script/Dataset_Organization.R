
setwd('C:\\Users/rli3/Documents/Publications/CancerMIRNome/')
library(RSQLite)

ccma.primary <- read.delim('shinyApp/data/miRNomes_Datasets_Primary.txt', header = T, sep = '\t', stringsAsFactors = F)
View(ccma.primary)

saveRDS(ccma.primary, file='shinyApp/data/miRNomes_Datasets_Primary.RDS')

ccma <- read.delim('shinyApp/data/miRNomes_Datasets.txt', header = T, sep = '\t', stringsAsFactors = F)
View(ccma)

saveRDS(ccma, file='shinyApp/data/miRNomes_Datasets.RDS')


meta.ccma <- list()
expr.ccma <- list()

for (i in 1:nrow(ccma)) {
  
  accession <- ccma$Accession[i]
  platform <- ccma$Annotation[i]
  dataset <- ccma$Dataset[i]
  
  print (i)
  print (dataset)
  
  exprData <- readRDS(file=paste0('shinyApp/data/Datasets/', accession, '_', platform, '_Expression.RDS'))
  phenoData <- readRDS(file=paste0('shinyApp/data/Datasets/', accession, '_', platform, '_Metadata.RDS'))

  expr.ccma[[dataset]] <- exprData
  meta.ccma[[dataset]] <- phenoData
  
  
}

saveRDS(expr.ccma, file='shinyApp/data/miRNomes_Expression.RDS')
saveRDS(meta.ccma, file='shinyApp/data/miRNomes_Metadata.RDS')


ccma$Dataset[1:2]

meta.ccma[[ccma$Dataset[1]]]

expr.ccma <- readRDS(file='shinyApp/data/miRNomes_Expression.RDS')
x <- lapply(expr.ccma, nrow)
x[which(x>1000 & x<2000)]
length(x[which(x>1000 & x<2000)])


x[which(x<1000)]
length(x[which(x<1000)])


names(expr.ccma)



library(RSQLite)
db <- dbConnect(SQLite(), dbname='CCMA_Expression.sqlite')
dbListTables(db)

for (project in names(expr.ccma)) {
  print (project)
  
  expr.table <- data.frame(expr.ccma[[project]], stringsAsFactors = F)
  #colnames(expr.table) <- gsub('.', '-', colnames(expr.table), fixed = T)
  #expr.table <- data.frame(expr.tcga.list[[project]][,1:900], stringsAsFactors = F)
  
  if (ncol(expr.table)>900) {
    saveRDS(expr.table, file=paste0('shinyApp/data/CCMA_Expression.', project, '.RDS'))
    next
    
  }
  
  dbWriteTable(db, project, expr.table, row.names = TRUE)
}


dbListTables(db)

dbDisconnect(db)


getCCMATable <- function(project) {
  db <- dbConnect(SQLite(), dbname='CCMA_Expression.sqlite')
  expr.table <- dbReadTable(db, project, row.names=TRUE)
  dbDisconnect(db)
  
  return (expr.table)
  
}

test <- getCCMATable('GSE39845')
View(test)

projects <- dbListTables(db)
projects

expr.data <- list()
for (prj in projects) {
  expr.data[[prj]] <- getCCMATable(prj)
  
}










#################

ccma <- readRDS(file='shinyApp/data/CCMA_Datasets.RDS')
idx

i <- 66
ccma$Accession[i]
ccma$Annotation[i]
expr <- readRDS(paste0('shinyApp/data/Datasets/',ccma$Accession[i], '_', 
               ccma$Annotation[i], '_Expression.RDS'))
genes <- rownames(expr)
genes

expr[1:5,1:5]

min(expr, na.rm = T)
max(expr, na.rm = T)


expr[expr<=0.1] <- 0.1

expr <- apply(expr, 2, as.numeric)

expr <- apply(expr, 2, function(v) round(v,3))
expr[1:5,1:5]

expr <- data.frame(expr, stringsAsFactors = F, row.names = genes)

expr <- apply(expr, 2, function(v) round(log2(v),3))
expr[1:5,1:5]

filter <- which(rowSums(is.na(expr))>0.2*ncol(expr))
filter
length(filter)
nrow(expr)

expr <- expr[-filter,]

saveRDS(expr, paste0('shinyApp/data/Datasets/',ccma$Accession[i], '_', 
               ccma$Annotation[i], '_Expression.RDS'))


idx <- c()

for (i in 1:nrow(ccma)) {
  ccma$Accession[i]
  ccma$Annotation[i]
  expr <- readRDS(paste0('shinyApp/data/Datasets/',ccma$Accession[i], '_', 
                         ccma$Annotation[i], '_Expression.RDS'))
  
  filter <- which(rowSums(is.na(expr))>0.2*ncol(expr))
  
  if (length(filter) > 0) {
    idx <- c(idx, i)
    print (i)
  }
  
}



i <- 38
ccma$Accession[i]
ccma$Annotation[i]
expr <- readRDS(paste0('shinyApp/data/Datasets/',ccma$Accession[i], '_', 
                       ccma$Annotation[i], '_Expression.RDS'))

filter <- which(colSums(is.na(expr))>0.5*nrow(expr))
filter

expr <- expr[,-filter]

pheno <- readRDS(paste0('shinyApp/data/Datasets/',ccma$Accession[i], '_', 
                        ccma$Annotation[i], '_Metadata.RDS'))
pheno

pheno <- pheno[-filter,]


saveRDS(pheno, paste0('shinyApp/data/Datasets/',ccma$Accession[i], '_', 
                     ccma$Annotation[i], '_Metadata.RDS'))

nrow(pheno)




min.val <- c()
max.val <- c()
na.val <- c()
char.val <- c()

for (i in 1:nrow(ccma)) {
  expr <- readRDS(paste0('shinyApp/data/Datasets/',ccma$Accession[i], '_', 
                         ccma$Annotation[i], '_Expression.RDS'))
  
  char.val <- c(char.val, sum(!is.numeric(expr)))
  expr <- apply(expr, 2, as.numeric)
  
  min.val <- c(min.val, min(expr, na.rm = T))
  max.val <- c(max.val, max(expr, na.rm = T))
  na.val <- c(na.val, sum(is.na(expr)))
  
}

which(max.val>100)
which(na.val>0)

which(char.val>0)





for (i in 1:nrow(ccma)) {
  message(i)
  ccma$Accession[i]
  ccma$Annotation[i]
  expr <- readRDS(paste0('shinyApp/data/Datasets/',ccma$Accession[i], '_', 
                         ccma$Annotation[i], '_Expression.RDS'))
  
  print (nrow(expr))
  idx <- which(rowSums(is.na(expr))>0)
  print (length(idx))
  
  
}






###########################################3


projects <- read.table('data/fromTCGA/projects-table.2020-06-20.tsv', header = T, sep = '\t', stringsAsFactors = F)
projects <- projects[1:33,]

exprLogCPM <- readRDS(file=paste0('data/fromTCGA/RNAseq_LogCPM_', projects$Project[2], '.All.RDS'))
genes <- intersect(gene.ids$Target.Ensembl, rownames(exprLogCPM))
genes
nrow(gene.ids)
length(genes)


expr.tcga.list <- list()

for (project in projects$Project[1:33]) {
  
  exprLogCPM <- readRDS(file=paste0('data/fromTCGA/RNAseq_LogCPM_', project, '.All.RDS'))
  
  exprLogCPM <- apply(exprLogCPM, 2, function(v) round(v,3))
  
  expr.tcga.list[[project]] <- exprLogCPM[gene.ids$Target.Ensembl,]
  
}

saveRDS(expr.tcga.list, file='data/fromTCGA/RNAseq_Expression_TCGA.miRTarBase.RDS')


##############
expr.tcga.list <- readRDS(file='shinyApp/data/RNAseq_Expression_TCGA.miRTarBase.RDS')
names(expr.tcga.list)

saveRDS(expr.tcga.list[['TCGA-BRCA']], file='RNAseq_Expression_TCGA.miRTarBase.TCGA-BRCA.RDS')


library(RSQLite)
db <- dbConnect(SQLite(), dbname='RNAseq_Expression_TCGA.miRTarBase.sqlite')
dbListTables(db)

for (project in names(expr.tcga.list)[2:33]) {
  print (project)
  
  expr.table <- data.frame(expr.tcga.list[[project]], stringsAsFactors = F)
  colnames(expr.table) <- gsub('.', '-', colnames(expr.table), fixed = T)
  #expr.table <- data.frame(expr.tcga.list[[project]][,1:900], stringsAsFactors = F)
  dbWriteTable(db, project, expr.table, row.names = TRUE)
}


dbListTables(db)
dbListFields(db, "TCGA-BLCA")

expr.table <- dbReadTable(db, "TCGA-BLCA", row.names=FALSE)
expr.table

datatable(as.data.frame(expr.table)[1:5,], 
          options = list(scrollX = TRUE, pageLength = 5))

dbDisconnect(db)


getRNATable <- function(project) {
  db <- dbConnect(SQLite(), dbname='RNAseq_Expression_TCGA.miRTarBase.sqlite')
  expr.table <- dbReadTable(db, project, row.names=TRUE)
  dbDisconnect(db)
  
  return (expr.table)
  
}

test <- getRNATable('TCGA-CHOL')
View(test)




#####
tcga <- read_excel('shinyApp/data/TCGA_Projects.xlsx')
tcga
tcga <- data.frame(tcga, stringsAsFactors = F)
tcga

saveRDS(tcga, file='shinyApp/data/TCGA_Projects.RDS')



###############################################

#### Additional clinical data

meta.tcga <- readRDS('shinyApp/data/Metadata_TCGA.RDS')

for (prj in names(meta.tcga)) {
  meta <- meta.tcga[[prj]]
  
  if (prj=='TCGA-PRAD') {
    
    prad <- readRDS('shinyApp/data/fromTCGA/Clinical_TCGA_PRAD_With_PreopPSA_and_BCR.RDS')
    
    samples <- intersect(rownames(prad), rownames(meta))
    meta$gleason_score <- meta$preop_psa <- NA
    
    meta[samples,]$gleason_score <- paste0(prad[samples,]$gleason_score, ' (', prad[samples,]$primary_pattern,
                                           '+', prad[samples,]$secondary_pattern, ')')
    
    meta[samples,]$preop_psa <- prad[samples,]$preop_psa
    
  } else {

    meta$preop_psa <- meta$gleason_score <- NA
    
  }
  
  meta[is.na(meta)] <- 'NA'
  meta.tcga[[prj]] <- meta
  
}


saveRDS(meta.tcga, 'shinyApp/data/Metadata_TCGA.RDS')

