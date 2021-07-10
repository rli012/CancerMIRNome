setwd('~/CancerMIRNome/')

mirtarbase <- readRDS(file='shinyApp/data/miRTarBase_hsa_MTI.RDS')
filter <- which(duplicated(paste(mirtarbase$miRBase.Accession, mirtarbase$Target.Ensembl)))
filter

# which(duplicated(paste(mirtarbase$miRBase.Accession, mirtarbase$Target.Entrez)[-filter]))

mirtarbase <- mirtarbase[-filter,]
dim(mirtarbase)

length(unique(mirtarbase$miRBase.Accession))


filter <- which(duplicated(mirtarbase$Target.Ensembl))
filter
gene.ids <- mirtarbase[-filter,]
gene.ids

which(!gene.ids$Target.Ensembl %in% rownames(rna.tcga[[1]]))


mirbase <- readRDS('shinyApp/data/miRBase_10.0_22.RDS')
mirbase

meta.tcga <- readRDS('shinyApp/data/Metadata_TCGA.RDS')
mir.tcga <- readRDS('shinyApp/data/miRNA_Expression_TCGA.RDS')
rna.tcga <- readRDS('shinyApp/data/RNAseq_Expression_TCGA.miRTarBase.RDS')

projects <- names(meta.tcga)
projects


cor.list <- list()
#cor.expr <- list()

for (prj in projects) {
  
  print (prj)
  
  cor.list[[prj]] <- list()
  
  rna.expr <- rna.tcga[[prj]]
  mir.expr <- mir.tcga[[prj]]
  
  samples <- intersect(colnames(rna.expr), colnames(mir.expr))
  
  rna.expr <- rna.expr[,samples]
  mir.expr <- mir.expr[,samples]
  
  #group <- meta.tcga[[prj]][samples,'sample_type']
  
  #cor.expr[[prj]] <- data.frame(t(mir.expr), t(rna.expr), Sample.Type=group,
  #                              Project=prj, stringsAsFactors = F)
  
  i <- 0
  for (mir in mirbase$ID) {

    i <- i + 1
    print (i)
    
    idx <- which(mirtarbase$miRBase.Accession==mir)
    
    if (length(idx)==0) {
      
      miRNA.Accession <- miRNA.ID <- Target.Ensembl <- 
        Target.Symbol <- Correlation <- P.Value <- BH.Adj.P <- character()
      
      cor.list[[prj]][[mir]] <- data.frame(miRNA.Accession, miRNA.ID, Target.Ensembl, 
                                Target.Symbol, Correlation, P.Value, BH.Adj.P,
                                stringsAsFactors = F)
      
      next
      
    }
    
    if (! mir %in% rownames(mir.expr)) {
      
      miRNA.Accession <- miRNA.ID <- Target.Ensembl <- 
        Target.Symbol <- Correlation <- P.Value <- BH.Adj.P <- character()
      
      cor.list[[prj]][[mir]] <- data.frame(miRNA.Accession, miRNA.ID, Target.Ensembl, 
                                    Target.Symbol, Correlation, P.Value, BH.Adj.P,
                                    stringsAsFactors = F)
      
      next
      
    }
    
    targets <- mirtarbase$Target.Ensembl[idx]
    
    #print (which(!targets %in% rownames(rna.expr)))
    targets <- intersect(targets, rownames(rna.expr))
    
    cor.val <- apply(rna.expr[targets,], 1, function(v) 
      cor.test(mir.expr[mir,], v, alternative='two.sided', method = 'pearson')$estimate)
    
    cor.p <- apply(rna.expr[targets,], 1, function(v) 
      cor.test(mir.expr[mir,], v, alternative='two.sided', method = 'pearson')$p.value)
    
    cor.fdr <- p.adjust(cor.p, method = 'BH')
    
    mir.name <- mirbase$Name[mirbase$ID==mir]
    targets.name <- gene.ids$Target.Symbol[match(targets, gene.ids$Target.Ensembl)]
    
    cor.data <- data.frame(mir, mir.name, targets, targets.name, cor.val, cor.p, cor.fdr, stringsAsFactors = F)

    #keep <- which(cor.val<0)
    
    #cor.data <- cor.data[keep,]
    
    o <- order(cor.data$cor.fdr, decreasing = F)
    
    cor.data <- cor.data[o,]
    
    cor.data$cor.val <- formatC(cor.data$cor.val, format = 'f', digits=3)
    cor.data$cor.p <- formatC(cor.data$cor.p, format = 'e', digits=2)
    cor.data$cor.fdr <- formatC(cor.data$cor.fdr, format = 'e', digits=2)
    
    colnames(cor.data) <- c('miRNA.Accession','miRNA.ID','Target.Ensembl','Target.Symbol','Correlation','P.Value','BH.Adj.P')
    rownames(cor.data) <- NULL
    
    cor.list[[prj]][[mir]] <- cor.data
    
  }

}

saveRDS(cor.list, file='shinyApp/data/Pearson.Correlation.miRTarBase.RDS')

cor.list <- readRDS(file='shinyApp/data/Pearson.Correlation.miRTarBase.RDS')
names(cor.list)
cor.list[[1]]

library(RSQLite)
db <- dbConnect(SQLite(), dbname='Pearson.Correlation.miRTarBase.sqlite')
dbListTables(db)
names(cor.list)

for (project in names(cor.list)) {
  print (project)
  
  cor.table <- do.call(rbind, cor.list[[project]])
  
  dbWriteTable(db, project, cor.table, row.names = FALSE)
}


dbListTables(db)
dbListFields(db, "TCGA-BLCA")

cor.table <- dbReadTable(db, "TCGA-BLCA", row.names=FALSE)
cor.table

datatable(as.data.frame(cor.table)[1:5,], 
          options = list(scrollX = TRUE, pageLength = 5))

dbDisconnect(db)



getCorTable <- function(project, mir) {
  db <- dbConnect(SQLite(), dbname='Correlation.miRTarBase.sqlite')
  cor.table <- dbReadTable(db, project, row.names=FALSE)
  dbDisconnect(db)
  
  return (cor.table)
  
}


test <- getCorTable('TCGA-BLCA')
test



getCorTable <- function(project, mir) {

  db <- dbConnect(SQLite(), dbname='Correlation.miRTarBase.sqlite')
  
  CMD <- paste0("SELECT * FROM [", project, "] WHERE [miRNA.Accession] == '", mir, "'")
  
  query <- dbSendQuery(db, CMD)
  cor.table <- dbFetch(query)
  
  dbClearResult(query)
  dbDisconnect(db)
  
  return (cor.table)
  
}


project <- 'TCGA-BLCA'
mir <- 'MIMAT0000062'
CMD <- paste0("SELECT * FROM [", project, "] WHERE [miRNA.Accession] == '", mir, "'")
CMD

res <- dbSendQuery(db, CMD)

res <- dbSendQuery(db, "SELECT * FROM [TCGA-BLCA] WHERE [miRNA.Accession] == 'MIMAT0000062'")
test <- dbFetch(res)

dbClearResult(res)

test

dbDisconnect(db)

test <- getCorTable(project = 'TCGA-BLCA', mir = 'MIMAT0000062')
test



########## SLOW

db <- dbConnect(SQLite(), dbname='shinyApp/data/Correlation.miRTarBase.sqlite')
dbListTables(db)

for (project in names(cor.list)) {
  print (project)
  
  for (mir in names(cor.list[[project]])) {
    cor.table <- cor.list[[project]][[mir]]
    dbWriteTable(db, paste0(project,'.',mir), cor.table, row.names = FALSE)
  }

}


dbListTables(db)
dbListFields(db, "TCGA-BLCA.MIMAT0015377")

cor.table <- dbReadTable(db, "TCGA-BLCA.MIMAT0015377", row.names=FALSE)
cor.table


datatable(as.data.frame(cor.table)[1:5,], 
          options = list(scrollX = TRUE, pageLength = 5))

dbDisconnect(db)


getCorTable <- function(project, mir) {
  db <- dbConnect(SQLite(), dbname='shinyApp/data/Correlation.miRTarBase.sqlite')
  cor.table <- dbReadTable(db, paste0(project,'.',mir), row.names=FALSE)
  
  dbDisconnect(db)
  
  return (cor.table)
  
  
  
}

test <- getCorTable('TCGA-BLCA', 'MIMAT0015377')
test




