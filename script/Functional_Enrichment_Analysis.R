
mirtarbase <- read_excel('shinyApp/data/Targets/mirTarBase_hsa_MTI.xlsx')
mirtarbase[1:5,1:5]
mirtarbase$miRNA[1:5]


mirbase <- readRDS('shinyApp/data/miRBase_10.0_22.RDS')
mirbase

entrez.ensembl.symbol <- readRDS('shinyApp/data/Homo_Sapiens_Gene_Annotation_ENSEMBL_HGNC_ENTREZ.RDS')
entrez.ensembl.symbol

idx <- match(mirtarbase$miRNA, mirbase$Name)
idx

mirtarbase$miRBase.Accession <- mirbase$ID[idx]
View(mirtarbase)


idx <- match(mirtarbase$`Target Gene (Entrez Gene ID)`, entrez.ensembl.symbol$entrez_id)
idx

mirtarbase$Ensembl <- entrez.ensembl.symbol$ensembl_id[idx]

mirtarbase <- data.frame(mirtarbase[,c(1:2,10,3:5,11,6:9)], stringsAsFactors = F)

colnames(mirtarbase) <- c('miRTarBase.ID','miRBase.ID','miRBase.Accession','miRNA.Species',
                          'Target.Symbol','Target.Entrez','Target.Ensembl','Target.Species',
                          'Experiments','Support.Type','References.PMID')


saveRDS(mirtarbase, file='shinyApp/data/miRTarBase_hsa_MTI.RDS')


####################################################################################

library(clusterProfiler)
library(GSEABase)
library(ReactomePA)
library(DOSE)

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




KEGG <- list()
i <- 0

for (mir in mirbase$ID) {
  i <- i + 1
  print (i)
  
  idx <- which(mirtarbase$miRBase.Accession==mir)
  
  if (length(idx)==0) {
    
    ID <- Description <- Count <-  List.Total <- Pop.Hits <- 
      Pop.Total <- Fold.Enrichment <- P.Value <- BH.Adj.P <- 
      Entrez <- Ensembl <- Symbol <- character()
    
    KEGG[[mir]] <- data.frame(ID,Description,
                              Count, List.Total, Pop.Hits, Pop.Total,
                              Fold.Enrichment, P.Value, BH.Adj.P,
                              Entrez, Ensembl, Symbol,
                              stringsAsFactors = F)
    
    next
  }
  
  targets <- mirtarbase$Target.Entrez[idx]
  
  
  
  # ### KEGG
  # kk <- enrichKEGG(gene          = targets,
  #                  organism      = 'hsa',
  #                  pAdjustMethod = 'BH',
  #                  pvalueCutoff  = 0.05,
  #                  minGSSize     = 10,
  #                  maxGSSize     = 500)
  # 
  # ### ReactomePA
  kk <- enrichPathway(gene          = targets,
                      organism      = 'human',
                      pAdjustMethod = 'BH',
                      pvalueCutoff  = 0.05,
                      minGSSize     = 10,
                      maxGSSize     = 500,
                      readable      = FALSE)
  
  
  ## Disease Ontology
  # kk <- enrichDO(gene          = targets,
  #                ont           = "DO",
  #                pAdjustMethod = 'BH',
  #                pvalueCutoff  = 0.05,
  #                minGSSize     = 10,
  #                maxGSSize     = 500,
  #                readable      = FALSE)
  # 
  # # ### Network of Cancer Gene
  # kk <- enrichNCG(gene        = targets,
  #                 pAdjustMethod = 'BH',
  #                 pvalueCutoff  = 0.05,
  #                 minGSSize     = 10,
  #                 maxGSSize     = 500,
  #                 readable      = FALSE)
  # 
  # ### DisGeNET
  # kk <- enrichDGN(gene          = targets,
  #                 pAdjustMethod = 'BH',
  #                 pvalueCutoff  = 0.05,
  #                 minGSSize     = 10,
  #                 maxGSSize     = 500,
  #                 readable      = FALSE)
  # 
  # ### Gene Ontology, BP
  # kk <- enrichGO(gene          = targets,
  #                OrgDb         = org.Hs.eg.db,
  #                ont           = "BP",
  #                pAdjustMethod = "BH",
  #                pvalueCutoff  = 0.05,
  #                minGSSize     = 10,
  #                maxGSSize     = 500,
  #                readable      = FALSE)
  #                 
  # ### Gene Ontology, CC
  # kk <- enrichGO(gene          = targets,
  #                OrgDb         = org.Hs.eg.db,
  #                ont           = "CC",
  #                pAdjustMethod = "BH",
  #                pvalueCutoff  = 0.05,
  #                minGSSize     = 10,
  #                maxGSSize     = 500,
  #                readable      = FALSE)
  # 
  # 
  # ### Gene Ontology, MF
  # kk <- enrichGO(gene          = targets,
  #                OrgDb         = org.Hs.eg.db,
  #                ont           = "MF",
  #                pAdjustMethod = "BH",
  #                pvalueCutoff  = 0.05,
  #                minGSSize     = 10,
  #                maxGSSize     = 500,
  #                readable      = FALSE)
  # 
  # 
  # ### MSigDb, H: hallmark gene sets
  # gmtfile <- 'shinyApp/data/MSigDB/h.all.v7.1.entrez.gmt'
  # h <- read.gmt(gmtfile)
  # 
  # kk <- enricher(gene          = targets,
  #                TERM2GENE     = h,
  #                pAdjustMethod = "BH",
  #                pvalueCutoff  = 0.05,
  #                minGSSize     = 10,
  #                maxGSSize     = 500)
  # 
  # ### MSigDb, C4: computational gene sets -> CGN: cancer gene neighborhoods
  # gmtfile <- 'shinyApp/data/MSigDB/c4.cgn.v7.1.entrez.gmt'
  # h <- read.gmt(gmtfile)
  # 
  # kk <- enricher(gene          = targets,
  #                TERM2GENE     = h,
  #                pAdjustMethod = "BH",
  #                pvalueCutoff  = 0.05,
  #                minGSSize     = 10,
  #                maxGSSize     = 500)
  # 
  # ### MSigDb, C4: computational gene sets -> CM: cancer modules
  # gmtfile <- 'shinyApp/data/MSigDB/c4.cm.v7.1.entrez.gmt'
  # h <- read.gmt(gmtfile)
  # 
  # kk <- enricher(gene          = targets,
  #                TERM2GENE     = h,
  #                pAdjustMethod = "BH",
  #                pvalueCutoff  = 0.05,
  #                minGSSize     = 10,
  #                maxGSSize     = 500)
  # 
  # ### MSigDb, C6: oncogenic signatures
  # gmtfile <- 'shinyApp/data/MSigDB/c6.all.v7.1.entrez.gmt'
  # h <- read.gmt(gmtfile)
  # 
  # kk <- enricher(gene          = targets,
  #                TERM2GENE     = h,
  #                pAdjustMethod = "BH",
  #                pvalueCutoff  = 0.05,
  #                minGSSize     = 10,
  #                maxGSSize     = 500)
  # 
  # ### MSigDb, C7: immunologic signatures
  # gmtfile <- 'shinyApp/data/MSigDB/c7.all.v7.1.entrez.gmt'
  # h <- read.gmt(gmtfile)
  # 
  # kk <- enricher(gene          = targets,
  #                TERM2GENE     = h,
  #                pAdjustMethod = "BH",
  #                pvalueCutoff  = 0.05,
  #                minGSSize     = 10,
  #                maxGSSize     = 500)
  
  
  
  if (is.null(kk)) {
    
    ID <- Description <- Count <-  List.Total <- Pop.Hits <- 
      Pop.Total <- Fold.Enrichment <- P.Value <- BH.Adj.P <- 
      Entrez <- Ensembl <- Symbol <- character()
    
    KEGG[[mir]] <- data.frame(ID,Description,
                              Count, List.Total, Pop.Hits, Pop.Total,
                              Fold.Enrichment, P.Value, BH.Adj.P,
                              Entrez, Ensembl, Symbol,
                              stringsAsFactors = F)
    next
  }
  
  if (nrow(kk)==0) {
    
    ID <- Description <- Count <-  List.Total <- Pop.Hits <- 
      Pop.Total <- Fold.Enrichment <- P.Value <- BH.Adj.P <- 
      Entrez <- Ensembl <- Symbol <- character()
    
    KEGG[[mir]] <- data.frame(ID,Description,
                              Count, List.Total, Pop.Hits, Pop.Total,
                              Fold.Enrichment, P.Value, BH.Adj.P,
                              Entrez, Ensembl, Symbol,
                              stringsAsFactors = F)
    next
  }
  
  List.Total <- length(targets)
  Count <- kk$Count
  Pop.Hits <- as.numeric(sapply(kk$BgRatio, function(x) strsplit(x, '/')[[1]][1]))
  Pop.Total <- as.numeric(sapply(kk$BgRatio, function(x) strsplit(x, '/')[[1]][2]))
  
  Fold.Enrichment <- format(Count/List.Total*Pop.Total/Pop.Hits, digits=2, nsmall=2)
  
  P.Value <- format(kk$pvalue, digits=3, nsmall=3)
  BH.Adj.P <- formatC(kk$p.adjust, format = 'e', digits=2)
  
  Entrez <- kk$geneID
  genes <- strsplit(kk$geneID, '/')
  
  idx <- lapply(genes, function(v) match(v, gene.ids$Target.Entrez))
  Ensembl <- unlist(lapply(idx, function(x) paste(gene.ids$Target.Ensembl[x], collapse = '/')))
  Symbol <- unlist(lapply(idx, function(x) paste(gene.ids$Target.Symbol[x], collapse = '/')))
  
  KEGG[[mir]] <- data.frame(ID = kk$ID,
                            Description = kk$Description,
                            Count, List.Total, Pop.Hits, Pop.Total,
                            Fold.Enrichment, P.Value, BH.Adj.P,
                            Entrez, Ensembl, Symbol,
                            stringsAsFactors = F)
  
}


saveRDS(KEGG, file='shinyApp/data/miRTarBase.KEGG.RDS')
saveRDS(KEGG, file='shinyApp/data/miRTarBase.REACTOME.RDS') #

saveRDS(KEGG, file='shinyApp/data/miRTarBase.DO.RDS')
saveRDS(KEGG, file='shinyApp/data/miRTarBase.NCG.RDS')
saveRDS(KEGG, file='shinyApp/data/miRTarBase.DGN.RDS') #

saveRDS(KEGG, file='shinyApp/data/miRTarBase.GOBP.RDS') #
saveRDS(KEGG, file='shinyApp/data/miRTarBase.GOCC.RDS') #
saveRDS(KEGG, file='shinyApp/data/miRTarBase.GOMF.RDS') #


saveRDS(KEGG, file='shinyApp/data/miRTarBase.MSigDBHALLMARK.RDS')
saveRDS(KEGG, file='shinyApp/data/miRTarBase.MSigDBC4CGN.RDS')
saveRDS(KEGG, file='shinyApp/data/miRTarBase.MSigDBC4CM.RDS')
saveRDS(KEGG, file='shinyApp/data/miRTarBase.MSigDBC6.RDS')
saveRDS(KEGG, file='shinyApp/data/miRTarBase.MSigDBC7.RDS') #

genesets <- c('KEGG','REACTOME','DO','NCG','DGN',
              'GOBP','GOCC','GOMF',
              'MSigDBHALLMARK','MSigDBC4CGN','MSigDBC4CM',
              'MSigDBC6','MSigDBC7')

enrichment.list <- list()
for (geneset in genesets) {
  
  enrichment.list[[geneset]] <- readRDS(paste0('shinyApp/data/miRTarBase.',geneset,'.RDS'))
  
  enrichment.list[[geneset]] <- lapply(enrichment.list[[geneset]],
                                       function(x) {
                                         x$Symbol <- gsub('/', '; ', x$Symbol)
                                         return(x)
                                       }
                                       )
    
  
  
} 

saveRDS(enrichment.list, file='shinyApp/data/Enrichment.miRTarBase.RDS')


View(enrichment.list[[geneset]])
enrichment.list[[geneset]][[9]]


enrichment.list <- readRDS(file='shinyApp/data/Enrichment.miRTarBase.RDS')
enrichment.list

names(enrichment.list)

names(enrichment.list[['KEGG']])
enrich.table <- enrichment.list[['KEGG']][['MIMAT0000062']]
enrich.table$Count <- paste0(enrich.table$Count, '/', enrich.table$List.Total)

enrich.table$Pop.Hits <- paste0(enrich.table$Pop.Hits, '/', enrich.table$Pop.Total)

colnames(enrich.table)[c(3,5)] <- c('Count/List.Total','Pop.Hits/Pop.Total')
enrich.table <- enrich.table[-c(4,6,7)]
enrich.table
