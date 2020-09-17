################################## Server #####################################

library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyjs)
#library(shinyWidgets)
library(shinycssloaders)
#library(readxl)
library(ggplot2)
library(Matrix)
library(stringr)
library(Biobase)
library(survival)
library(survminer)
library(limma)
library(edgeR)
library(DT)
library(dplyr)
library(pROC)
library(ROCR)
library(digest)
library(pheatmap)
library(plotly)
library(glmnet)


library(dashboardthemes)
library(shinythemes)

source('script/shiny_functions.R')


################################## Data #####################################

### TCGA Datasets
tcga.datasets <- readRDS('data/TCGA_Projects.RDS')

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
ccma.primary <- readRDS('data/miRNomes_Datasets_Primary.RDS')


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


################################## Input #####################################

###### miRNA

mir.default <- 'MIMAT0000062' # hsa-let-7a-5p
mir.annotation <- readRDS('data/miRBase_10.0_22.RDS')

mir.id <- selectizeInput(inputId = "mir.id", label=h4(strong('Search a miRNA')), #list(h4('Search a miRNA:'), icon('search', 'fa-1.5x')),# h4(strong('miRNA'))
                         choices = NULL, selected = mir.default, #mir.default, 
                         multiple = FALSE, width = 400,
                         options = list(placeholder = 'e.g. hsa-miR-7a-5p', #  or MIMAT0000062
                                        server = TRUE, selectOnTab=TRUE,
                                        searchField = c('Name', 'ID', 'Previous_ID'),
                                        labelField = "Name",
                                        valueField = "ID",
                                        #maxOptions = 5,
                                        render = I("{option: function(item, escape) 
                                                   {var gene = '<div>' + '<strong>' + escape(item.Name) + '</strong>' + '<ul>';
                                                   gene = gene + '<li>' + 'Previous IDs: ' + item.Previous_ID + '</li>';
                                                   gene = gene + '<li>' + 'Accession: ' + item.ID + '</li>' + '</ul>' + '</div>';
                                                   return gene
                                                   }
                                                   }")
                         ))


###### TCGA Project

project.default <- 'TCGA-BLCA'

project.id <- selectizeInput(inputId = "project.id", label=h5(strong('TCGA Project:')),# h4(strong('miRNA'))
                             choices = NULL, selected = project.default, 
                             multiple = FALSE, width = 150,
                             options = list(placeholder = 'Select a project',
                                            server = TRUE, selectOnTab=TRUE
                             ))

project.id.cor <- selectizeInput(inputId = "project.id.cor", label=h5(strong('TCGA Project:')),# h4(strong('miRNA'))
                                 choices = NULL, selected = project.default, 
                                 multiple = FALSE, width = 150,
                                 options = list(placeholder = 'Select a project',
                                                server = TRUE, selectOnTab=TRUE
                                 ))

###### Gene Sets

gene.sets <- c('Kyoto Encyclopedia of Genes and Genomes (KEGG)' = 'KEGG',
               'REACTOME' = 'REACTOME',
               'Disease Ontology (DO)' = 'DO',
               'Network of Cancer Gene (NCG)' = 'NCG',
               'DisGeNET' = 'DGN', 
               'Gene Ontology - Biological Process (GO-BP)' = 'GOBP',
               'Gene Ontology - Cellular Component (GO-CC)' = 'GOCC',
               'Gene Ontology - Molecular Function (GO-MF)' = 'GOMF',
               'MSigDB - H:HALLMARK' = 'MSigDBHALLMARK',
               'MSigDB - C4:Cancer Gene Neighborhoods' = 'MSigDBC4CGN',
               'MSigDB - C4:Cancer Modules' = 'MSigDBC4CM',
               'MSigDB - C6:Oncogenic Signatures' = 'MSigDBC6',
               'MSigDB - C7:Immunologic Signatures' = 'MSigDBC7')

geneset.id.default <- gene.sets[1]
geneset.id <- selectizeInput(inputId = "geneset.id", label=h5(strong('Gene Sets:')),# h4(strong('miRNA'))
                             choices = NULL, selected = geneset.id.default, 
                             multiple = FALSE, width = 410,
                             options = list(placeholder = 'Select a gene set',
                                            server = TRUE, selectOnTab=TRUE
                             ))

################################## Server #####################################

server <- function(input, output, session) { 
  
  updateSelectizeInput(session, 'mir.id', choices = mir.annotation, selected = mir.default, server = TRUE)
  updateSelectizeInput(session, 'project.id', choices = projects.tcga, selected = project.default, server = TRUE)
  updateSelectizeInput(session, 'project.id.cor', choices = projects.tcga, selected = project.default, server = TRUE)
  updateSelectizeInput(session, 'geneset.id', choices = gene.sets, selected = geneset.id.default, server = TRUE)
  
  seed <- reactiveVal()
  
  # observe({
  #   if (input$navbar == "home") {
  #     js$refresh();
  #   }
  # })
  
  
  ################################################################
  ######################## Information ###########################
  
  observeEvent(input$mir.id, {
    
    output$mir.name <- renderUI({ 
      mir.id <- input$mir.id
      mir.name <- mir.annotation[mir.id, 'Name']
      mir.url <- paste0('http://www.mirbase.org/cgi-bin/mature.pl?mature_acc=', mir.id)
      mir.url <- a(mir.name, href = mir.url, target="_blank", style = "font-size: 150%;")
      tagList(mir.url)
    })
    
    output$mir.preid <- renderText({ 
      mir.id <- input$mir.id
      mir.preid <- mir.annotation[mir.id, 'Previous_ID']
      mir.preid <- paste0('Previous IDs: ', mir.preid)
      mir.preid
    })
    
    
    
    output$mir.info <- renderText({ 
      mir.id <- input$mir.id
      mir.name <- mir.annotation[mir.id, 'Name']
      mir.info <- paste0('Accession: ', mir.id)
      mir.info
    })
    
    output$mir.seq <- renderText({ 
      mir.id <- input$mir.id
      mir.seq <- mir.annotation[mir.id, 'Sequence']
      mir.seq <- paste0('Sequence: ', mir.seq)
      mir.seq
    })
    
    
    output$mir.targets <- renderUI({ 
      mir.id <- input$mir.id
      mir.name <- mir.annotation[mir.id, 'Name']
      mir.encori <- paste0('http://starbase.sysu.edu.cn/agoClipRNA.php?source=mRNA&flag=miRNA&clade=mammal&genome=human&assembly=hg19&miRNA=',
                           mir.name, '&clipNum=&deNum=&panNum=&proNum=&program=&target=')
      mir.encori <- a('ENCORI', href = mir.encori, target="_blank", style = "font-size: 100%;")
      
      mir.mirdb <- paste0('http://mirdb.org/cgi-bin/search.cgi?searchType=miRNA&full=mirbase&searchBox=',mir.id)
      mir.mirdb <- a('miRDB', href = mir.mirdb, target="_blank", style = "font-size: 100%;")
      
      mir.mirtarbase <- paste0('http://mirtarbase.cuhk.edu.cn/php/search.php?opt=b_mirna&org=hsa&bname=',mir.name)
      mir.mirtarbase <- a('miRTarBase', href = mir.mirtarbase, target="_blank", style = "font-size: 100%;")
      
      mir.targetscan <- paste0('http://www.targetscan.org/cgi-bin/targetscan/vert_72/targetscan.cgi?mirg=',mir.name)
      mir.targetscan <- a('TargetScan', href = mir.targetscan, target="_blank", style = "font-size: 100%;")
      
      mir.dianatarbase <- paste0('http://diana.imis.athena-innovation.gr/DianaTools/index.php?r=tarbase/index&mirnas=',mir.name)
      mir.dianatarbase <- a('Diana-TarBase', href = mir.dianatarbase, target="_blank", style = "font-size: 100%;")
      
      tagList("Targets:", mir.encori, mir.mirdb, mir.mirtarbase, mir.targetscan, mir.dianatarbase)
    })
    
    
    
  })
  
  
  
  #########################################################
  ######################## TCGA ###########################
  
  observeEvent(input$mir.id, {
    
    
    tcga.overview <- reactiveValues()
    
    
    output$tcga_boxplot <- renderPlot({
      
      mir.id <- input$mir.id
      mir.name <- mir.annotation[mir.id, 'Name']
      
      group <- unlist(lapply(meta.tcga, function(x) x[,'sample_type']))
      expr <- unlist(lapply(mir.tcga, function(x) x[mir.id,]))
      project <- unlist(lapply(meta.tcga, function(x) x[,'project_id']))
      
      dataForBoxPlot <- data.frame(expr, group, project, mir=mir.name)
      
      tcga.overview$box.data <- dataForBoxPlot
      
      p <- tcgaboxplotFun(dataForBoxPlot)
      
      tcga.overview$box.plot <- p
      
      p
    })
    
    
    output$tcga.box.summ.downbttn.csv <- downloadHandler(
      filename = function(){paste('box.csv', sep = '')},
      
      content = function(file){
        write.csv(tcga.overview$box.data, file, row.names = FALSE, quote = F)
      })
    
    output$tcga.box.summ.downbttn.png <- downloadHandler(
      filename = function(){paste('box.png', sep = '')},
      
      content = function(file){
        png(file, width = 1000, height = 500)
        print(tcga.overview$box.plot)
        dev.off()
      })
    
    output$tcga.box.summ.downbttn.pdf <- downloadHandler(
      filename = function(){paste('box.pdf', sep = '')},
      
      content = function(file){
        pdf(file, width = 10, height = 5)
        print(tcga.overview$box.plot)
        dev.off()
      })
    

    output$tcga_rocplot_forest <- renderPlot({
      
      mir.id <- input$mir.id
      mir.name <- mir.annotation[mir.id, 'Name']
      
      group <- unlist(lapply(meta.tcga, function(x) x[,'sample_type']))
      expr <- unlist(lapply(mir.tcga, function(x) x[mir.id,]))
      project <- unlist(lapply(meta.tcga, function(x) x[,'project_id']))
      
      
      dataForForestPlot <- c()
      
      for (prj in projects.tcga.sub) {
        
        idx <- which(project==prj)
        
        roc.test <- roc(group[idx], expr[idx], plot=FALSE, ci=TRUE, auc=TRUE)
        ci.auc <- roc.test$ci
        
        auc <- ci.auc[2]
        auc.ci.lower95 <- ci.auc[1]
        auc.ci.upper95 <- ci.auc[3]
        
        auc <- format(auc, digits = 2, nsmall=2)
        auc.ci.lower95 <- format(auc.ci.lower95, digits = 2, nsmall=2)
        auc.ci.upper95 <- format(auc.ci.upper95, digits = 2, nsmall=2)
        
        
        # pred <- prediction(expr[idx], group[idx])
        # perf <- performance(pred, measure = "tpr", x.measure = "fpr")
        # 
        # FPR <- perf@x.values[[1]]
        # TPR <- perf@y.values[[1]]
        # 
        # #df <- data.frame(FPR,TPR)
        # 
        # auc.test <- wilcox.test(FPR, TPR, alternative = 'two.sided')
        # pvalue <- formatC(auc.test$p.value, format = 'e', digits = 2)
        
        dataForForestPlot <- rbind(dataForForestPlot, 
                                   c(auc, auc.ci.lower95, auc.ci.upper95)) #, pvalue
        
      }
      
      dataForForestPlot <- apply(dataForForestPlot, 2, as.numeric)
      
      dataForForestPlot <- data.frame(dataForForestPlot,
                                      row.names = projects.tcga.sub,
                                      stringsAsFactors = F)
      
      colnames(dataForForestPlot) <- c('AUC','Lower95','Upper95') #,'P.Value'
      
      dataForForestPlot$mir <- mir.name
      
      dataForForestPlot$Project <- projects.tcga.sub
      o <- order(dataForForestPlot$AUC, decreasing = F)
      
      dataForForestPlot$Project <- factor(dataForForestPlot$Project,
                                          levels = dataForForestPlot$Project[o])
      
      tcga.overview$roc.forest.data <- dataForForestPlot[rev(o),]
      
      p <- tcgaROCForestplotFun(dataForForestPlot)
      
      tcga.overview$roc.forest.plot <- p
      
      p
    })
    
    
    output$tcga.roc.forest.downbttn.csv <- downloadHandler(
      filename = function(){paste('roc.forest.csv', sep = '')},
      
      content = function(file){
        write.csv(tcga.overview$roc.forest.data, file, row.names = FALSE, quote = F)
      })
    
    output$tcga.roc.forest.downbttn.png <- downloadHandler(
      filename = function(){paste('roc.forest.png', sep = '')},
      
      content = function(file){
        png(file, width = 1000, height = 500)
        print(tcga.overview$roc.forest.plot)
        dev.off()
      })
    
    output$tcga.roc.forest.downbttn.pdf <- downloadHandler(
      filename = function(){paste('roc.forest.pdf', sep = '')},
      
      content = function(file){
        pdf(file, width = 10, height = 5)
        print(tcga.overview$roc.forest.plot)
        dev.off()
      })
    

    output$tcga_km_forest <- renderPlot({
      
      mir.id <- input$mir.id
      mir.name <- mir.annotation[mir.id, 'Name']
      
      group <- unlist(lapply(meta.tcga, function(x) x[,'sample_type']))
      expr <- unlist(lapply(mir.tcga, function(x) x[mir.id,]))
      project <- unlist(lapply(meta.tcga, function(x) x[,'project_id']))
      
      os.time <- unlist(lapply(meta.tcga, function(x) as.numeric(x[,'OS.time'])/30))
      os.status <- unlist(lapply(meta.tcga, function(x) as.numeric(x[,'OS'])))
      
      dataForForestPlot <- c()
      
      for (prj in projects.tcga) {
        
        idx <- which(project==prj & group=='Tumor')
        
        surv.data <- data.frame(expr=expr[idx],
                                os.time=os.time[idx],
                                os.status=os.status[idx],
                                stringsAsFactors = F)
        
        km <- kmTest(exprDa=surv.data$expr, daysToDeath=surv.data$os.time, vitalStatus=surv.data$os.status)
        dataForForestPlot <- rbind(dataForForestPlot, as.numeric(format(km, digits=3)))
        
      }
      
      dataForForestPlot <- apply(dataForForestPlot, 2, as.numeric)
      
      dataForForestPlot <- data.frame(dataForForestPlot,
                                      row.names = projects.tcga,
                                      stringsAsFactors = F)
      
      colnames(dataForForestPlot) <- c('HR','Lower95','Upper95','P.Value')
      
      dataForForestPlot$mir <- mir.name
      
      dataForForestPlot$Project <- projects.tcga
      o <- order(dataForForestPlot$HR, decreasing = F)
      
      dataForForestPlot$Project <- factor(dataForForestPlot$Project,
                                          levels = dataForForestPlot$Project[o])
      
      tcga.overview$km.forest.data <- dataForForestPlot[rev(o),]
      
      p <- tcgaKMForestplotFun(dataForForestPlot)
      
      tcga.overview$km.forest.plot <- p
      
      p
      
    })
    
    
    output$tcga.km.forest.downbttn.csv <- downloadHandler(
      filename = function(){paste('km.forest.csv', sep = '')},
      
      content = function(file){
        write.csv(tcga.overview$km.forest.data, file, row.names = FALSE, quote = F)
      })
    
    output$tcga.km.forest.downbttn.png <- downloadHandler(
      filename = function(){paste('km.forest.png', sep = '')},
      
      content = function(file){
        png(file, width = 1000, height = 650)
        print(tcga.overview$km.forest.plot)
        dev.off()
      })
    
    output$tcga.km.forest.downbttn.pdf <- downloadHandler(
      filename = function(){paste('km.forest.pdf', sep = '')},
      
      content = function(file){
        pdf(file, width = 10, height = 6.5)
        print(tcga.overview$km.forest.plot)
        dev.off()
      })
    
    
    
    observeEvent(input$project.id, {
      
      output$tcga_violinplot <- renderPlot({
        
        mir.id <- input$mir.id
        mir.name <- mir.annotation[mir.id, 'Name']
        
        project <- input$project.id
        
        group <- meta.tcga[[project]][,'sample_type']
        expr <- mir.tcga[[project]][mir.id,]
        
        dataForBoxPlot <- data.frame(expr, group, mir.id, project,
                                     stringsAsFactors = F)
        
        p <- BoxPlotFun(dataForBoxPlot)
        p
        
      })
      
      output$tcga_rocplot <- renderPlot({
        
        mir.id <- input$mir.id
        mir.name <- mir.annotation[mir.id, 'Name']
        
        project <- input$project.id
        
        group <- meta.tcga[[project]][,'sample_type']
        expr <- mir.tcga[[project]][mir.id,]
        
        dataForROCPlot <- data.frame(expr, group)
        dataForROCPlot$group <- ifelse(dataForROCPlot$group=='Normal',0,1)
        
        p <- rocplotFun(dataForROCPlot)
        p
      })
      
      output$tcga_km_plot <- renderPlot({
        
        mir.id <- input$mir.id
        mir.name <- mir.annotation[mir.id, 'Name']
        
        project <- input$project.id
        
        group <- meta.tcga[[project]][,'sample_type']
        expr <- mir.tcga[[project]][mir.id,]
        
        os.time <- as.numeric(meta.tcga[[project]][,'OS.time'])/30
        os.status <- as.numeric(meta.tcga[[project]][,'OS'])
        
        dataForKMPlot <- data.frame(expr, os.time, os.status, mir.id, project,
                                    stringsAsFactors = F)
        
        p <- KMPlotFun(dataForKMPlot)
        p
        
      })
      
      
    })
    
    
    output$correlation <- DT::renderDataTable({
      
      # if (!exists('cor.table')) {
      #   cor.table <- readRDS('data/Correlation.miRTarBase.RDS')
      # }
      
      mir <- input$mir.id
      project <- input$project.id.cor
      
      cor.table <- getCorTable(project = project, mir = mir)
      cor.table
      
    }, 
    callback=JS('$("button.buttons-copy").css("background","lightskyblue");
                 $("button.buttons-copy").css("border","white");
                 $("button.buttons-csv").css("background","lightskyblue");
                 $("button.buttons-csv").css("border","white");
                 $("button.buttons-excel").css("background","lightskyblue");
                 $("button.buttons-excel").css("border","white");
                 return table;'),
    options = list(pageLength = 5, dom = 'lfrtBp', buttons = c('copy', 'csv', 'excel')), #i
    extensions = "Buttons",
    selection = list(mode='single', selected=1),
    server = FALSE
    )
    
    observeEvent(input$correlation_rows_selected, {
      output$cor_plot <- renderPlot({
        
        mir <- input$mir.id
        project <- input$project.id.cor
        
        cor.table <- getCorTable(project = project, mir = mir)
        
        idx <- input$correlation_rows_selected
        mir.id <- cor.table[idx, 'miRNA.Accession']
        mir.name <- cor.table[idx, 'miRNA.ID']
        
        target.id <- cor.table[idx, 'Target.Ensembl']
        target.name <- cor.table[idx, 'Target.Symbol']
        
        rna.tcga <- getRNATable(project)
        colnames(rna.tcga) <- gsub('.', '-', colnames(rna.tcga), fixed = T)
        
        samples <- intersect(colnames(mir.tcga[[project]]), colnames(rna.tcga))
        
        mir.expr <- as.numeric(mir.tcga[[project]][mir.id,samples])
        rna.expr <- as.numeric(rna.tcga[target.id,samples])
        
        group <- meta.tcga[[project]][samples,'sample_type']

        coef <- cor.table[idx, 'Correlation']
        p.val <- cor.table[idx, 'P.Value']
        
        dataForCorrPlot <- data.frame(mir.expr, rna.expr, group, project,
                                      mir.id, target.id, mir.name, target.name,
                                      coef, p.val, stringsAsFactors = F)
        
        p <- ExprCorrPlotFun(dataForCorrPlot)
        p
        
      })
      
    })
    
    
    
    output$enrichment <- DT::renderDataTable({
      
      mir <- input$mir.id
      geneset <- input$geneset.id
      
      enrich.table <- enrichment.table[[geneset]][[mir]]
      
      if (nrow(enrich.table)==0) {
        shinyjs::hide('enrichment_bar_plot')
        shinyjs::hide('enrichment_bubble_plot')
      } else {
        shinyjs::show('enrichment_bar_plot')
        shinyjs::show('enrichment_bubble_plot')
      }
      
      if (nrow(enrich.table)>0) {
        enrich.table$Count <- paste0(enrich.table$Count, '/', enrich.table$List.Total)
        enrich.table$Pop.Hits <- paste0(enrich.table$Pop.Hits, '/', enrich.table$Pop.Total)
      }
      
      enrich.table <- enrich.table[-c(4,6,7,10,11)]
      colnames(enrich.table)[c(3,4)] <- c('Count/List.Total','Pop.Hits/Pop.Total')
      
      enrich.table
      
    }, 
    callback=JS('$("button.buttons-copy").css("background","lightskyblue");
                 $("button.buttons-copy").css("border","white");
                $("button.buttons-csv").css("background","lightskyblue");
                $("button.buttons-csv").css("border","white");
                $("button.buttons-excel").css("background","lightskyblue");
                $("button.buttons-excel").css("border","white");
                return table;'),
    options = list(pageLength = 5, dom = 'lfrtBp', buttons = c('copy', 'csv', 'excel')), #i
    extensions = "Buttons",
    selection = list(mode='none', selected=1), ### === not selectable
    server = FALSE
    )
    
    
    output$enrichment_bar_plot <- renderPlot({
      
      mir <- input$mir.id
      geneset <- input$geneset.id
      
      req(nrow(enrichment.table[[geneset]][[mir]])>0)
      
      dataForBarPlot <- enrichment.table[[geneset]][[mir]]
      dataForBarPlot$BH.Adj.P <- as.numeric(dataForBarPlot$BH.Adj.P)
      dataForBarPlot$Count <- as.numeric(dataForBarPlot$Count)
      dataForBarPlot$Fold.Enrichment <- as.numeric(dataForBarPlot$Fold.Enrichment)
      
      if (nrow(dataForBarPlot)>30) {
        dataForBarPlot <- dataForBarPlot[1:30,]
      }
      
      p <- EnrichmentBarPlotFun(dataForBarPlot)
      p
      
    })
    
    output$enrichment_bubble_plot <- renderPlot({
      
      mir <- input$mir.id
      geneset <- input$geneset.id
      
      req(nrow(enrichment.table[[geneset]][[mir]])>0)
      
      dataForBubblePlot <- enrichment.table[[geneset]][[mir]]
      dataForBubblePlot$BH.Adj.P <- as.numeric(dataForBubblePlot$BH.Adj.P)
      dataForBubblePlot$Count <- as.numeric(dataForBubblePlot$Count)
      dataForBubblePlot$Fold.Enrichment <- as.numeric(dataForBubblePlot$Fold.Enrichment)
      
      if (nrow(dataForBubblePlot)>30) {
        dataForBubblePlot <- dataForBubblePlot[1:30,]
      }
      
      p <- EnrichmentBubblePlotFun(dataForBubblePlot)
      p
      
    })
    
  })
  
  
  ########################################################################
  ######################## GENE-LEVEL ANALYSIS ###########################
  
  output$browser_datasets <- DT::renderDataTable({ccma.primary},
                                                 options = list(pageLength = 6), #8
                                                 selection = list(mode='multiple', selected=c(3,6)) #,7,8
  )
  
  observeEvent(input$mir.id, {
    req(mir.id)
    observeEvent(input$browser_datasets_rows_selected, {
      
      plot_data <- reactive({
        
        mir.id <- input$mir.id
        
        idx <- sort(input$browser_datasets_rows_selected)
        datasets <- as.character(ccma.datasets[idx,'Dataset'])
        
        expr.ccma <- list()
        
        for (dt in datasets) {
          expr.ccma[[dt]] <- getCCMATable(dt)
        }
        
        lapply(datasets, function(dataset) 
        {group <- meta.ccma[[dataset]][,'Disease.Status']
        expr <- expr.ccma[[dataset]][mir.id,]
        expr <- as.numeric(expr)
        
        dataForViolinPlot <- data.frame(expr, group, dataset,
                                        stringsAsFactors = F)
        
        expr.med <- dataForViolinPlot %>% dplyr::group_by(group) %>% 
          dplyr::summarise(med=median(expr, na.rm=T), N=length(expr))
        
        o <- order(expr.med$med, decreasing = T)
        expr.med <- expr.med[o,]
        
        idx <- match(dataForViolinPlot$group, expr.med$group)
        n <- expr.med$N[idx]
        
        expr.med$group <- paste0(expr.med$group, ' (N=', expr.med$N, ')')
        dataForViolinPlot$group <- paste0(dataForViolinPlot$group, ' (N=', n, ')')
        
        if (sum(grepl('Healthy', expr.med$group))==1) {
          idx <- grep('Healthy',expr.med$group)
          group.levels <- c(expr.med$group[idx],expr.med$group[-idx])
        } else {
          group.levels <- expr.med$group
        }
        
        dataForViolinPlot$group <- factor(dataForViolinPlot$group, levels=group.levels)
        
        return (dataForViolinPlot)
        }
        )
        
      })
      
      output$multi_plot_ui <- renderUI({
        
        idx <- sort(input$browser_datasets_rows_selected)
        datasets <- as.character(ccma.datasets[idx,'Dataset'])
        
        plot.widths <- lapply(datasets, function(dataset) 
        {group <- meta.ccma[[dataset]][,'Disease.Status']
        
        if (length(unique(group))<=5) {
          plot.width <- reactive(500)
        } else if (length(unique(group))>5 & length(unique(group))<10){
          plot.width <- reactive(800)
        } else if (length(unique(group))>=10){
          plot.width <- reactive(1000)
        } else {
          plot.width <- reactive(100 * length(unique(group)))
        }
        return (plot.width)
        
        })
        
        
        lapply(seq_along(input$browser_datasets_rows_selected),
               function(n) {
                 return(plot_overlay_ui(paste0("n", n), height=500, width=plot.widths[[n]]()))
               })
      })
      
      lapply(seq_along(input$browser_datasets_rows_selected),
             function(i){
               callModule(plot_overlay_server,
                          paste0("n", i),
                          dataForViolinPlot = plot_data()[[i]])
             }
      )
      
      
    })
    
    
  })
  
  
  
  
  ###########################################################################
  ######################## DATASET-LEVEL ANALYSIS ###########################
  
  
  #################### TCGA
  
  output$tcga_datasets <- DT::renderDataTable({tcga.datasets},
                                              options = list(pageLength = 5),
                                              selection = list(mode='single', selected=2)
  )
  
  observeEvent(input$tcga_datasets_rows_selected, {
    
    idx <- input$tcga_datasets_rows_selected
    req(idx)
    
    project <- as.character(tcga.datasets[idx,'Project'])
    
    meta <- meta.tcga[[project]]
    expr <- mir.tcga[[project]]
    
    groups <- meta$sample_type
    group.levels <- c('Normal','Tumor') ###
    groups <- factor(groups, levels = group.levels)
    
    groups <- as.data.frame(table(groups), stringsAsFactors=F)
    
    group.names <- sapply(groups[,1], function(x) digest(x,algo='murmur32',seed=sample(1e9,1)))
    
    groups <- cbind(rep(NA, nrow(groups)), rep(NA, nrow(groups)), groups)
    
    groups[,1] <- sprintf(
      '<input type="radio" name="%s" value="%s"/>',
      group.names, 1)
    
    groups[,2] <- sprintf(
      '<input type="radio" name="%s" value="%s"/>',
      group.names, 2)
    
    colnames(groups) <- c('Control','Case','Groups','N')
    
    output$groups.tcga <- DT::renderDataTable(groups, rownames = FALSE, escape = FALSE, selection = 'none', server = FALSE,
                                              options=list(dom = 'tp', paging = TRUE, pageLength = 100, #ordering = FALSE,
                                                           initComplete = JS("
                                                                             function(setting, json) {
                                                                             $(this.api().table().container())
                                                                             .find('div.dataTables_paginate')
                                                                             .css('display', this.api().page.info().pages <= 1 ? 'none' : 'block');
                                                                             }"),
        drawCallback = JS("
                          function(settings) {
                          Shiny.unbindAll(this.api().table().node());
                          Shiny.bindAll(this.api().table().node());
                          }")
        ),
        callback = JS("
                      table.rows().every(function(i, tab, row) {
                      var $this = $(this.node());
                      //$(\"input[name='\" + this.data()[2] + \"']\").prop('checked', false);
                      //console.log($this.children()[0]);
                      //console.log($.parseHTML(this.data()[0])[0].name);
                      $this.attr('id', $.parseHTML(this.data()[0])[0].name); //one time hash value
                      //$this.attr('id', this.data()[2]); //Group Name
                      $this.addClass('shiny-input-radiogroup');
                      //console.log($this.prop('checked'));
                      });
                      Shiny.unbindAll(table.table().node());
                      Shiny.bindAll(table.table().node());
                      "
                      )
                      )
    
    
    
    
    
    output$sel = renderPrint({
      unlist(sapply(group.names, function(i) input[[i]]))
    })
    
    
    output$dataset_summary_tcga <- renderText({ 
      dataset_summary_tcga <- as.character(paste0(tcga.datasets[idx,'Project'], ': ', tcga.datasets[idx,'Study.Name']))
      dataset_summary_tcga
    })
    
    output$pie_sample_type_tcga <- renderPlotly({
      sample.freq <- table(meta$sample_type)
      dataForPiePlot <- data.frame(num=as.numeric(sample.freq), sam=names(sample.freq))
      
      p <- piePlotlyFun(dataForPiePlot)
      p
    })
    
    output$pie_pstage_tcga <- renderPlotly({
      keep <- which(meta$sample_type=='Tumor')
      sample.freq <- table(meta$pathologic_stage[keep])
      dataForPiePlot <- data.frame(num=as.numeric(sample.freq), sam=names(sample.freq))
      
      p <- piePlotlyFun(dataForPiePlot)
      p
    })
    
    output$pie_cstage_tcga <- renderPlotly({
      keep <- which(meta$sample_type=='Tumor')
      sample.freq <- table(meta$clinical_stage[keep])
      dataForPiePlot <- data.frame(num=as.numeric(sample.freq), sam=names(sample.freq))
      
      p <- piePlotlyFun(dataForPiePlot)
      p
    })
    
    output$histogram_age_tcga <- renderPlotly({
      keep <- which(meta$sample_type=='Tumor')
      dataForHistogram <- meta[keep,]
      dataForHistogram$x <- as.integer(dataForHistogram$age_at_initial_pathologic_diagnosis)
      
      p <- histPlotlyFun(dataForHistogram)
      p
      #p <- histogramFun(dataForHistogram)
      #ggplotly(p, height = 400, width = 350)
      
    })
    
    
    output$pie_os_status_tcga <- renderPlotly({
      
      keep <- which(meta$sample_type=='Tumor')
      sample.freq <- table(as.numeric(meta$OS[keep]))
      dataForPiePlot <- data.frame(num=as.numeric(sample.freq), sam=names(sample.freq))
      dataForPiePlot$sam <- ifelse(dataForPiePlot$sam==0, 'Alive', 'Dead')
      
      p <- piePlotlyFun(dataForPiePlot)
      p
    })
    
    
    output$km_os_time_tcga <- renderPlotly({
      
      idx <- input$tcga_datasets_rows_selected
      project <- as.character(tcga.datasets[idx,'Project'])
      
      keep <- which(meta$sample_type=='Tumor')
      
      
      daysToDeath <- as.numeric(meta$OS.time[keep])/30
      vitalStatus <- as.numeric(meta$OS[keep])
      
      dataForKMPlot <- data.frame(daysToDeath, vitalStatus)
      
      fit <- survfit(Surv(daysToDeath, vitalStatus) ~ 1, data=dataForKMPlot)
      
      p <- ggsurvplot(fit, data=dataForKMPlot, #pval = paste(label1, '\n', label2), pval.coord = c(xpos, ypos1),
                      #pval.size=4,
                      font.main = c(12, 'bold', 'black'), conf.int = FALSE, 
                      #title = project,
                      legend = 'none', 
                      #color = c('blue', 'green'),
                      palette= c(google.blue, google.red),
                      #legend.labs = c(paste('Low Expr (N=',nL,')',sep=''), 
                      #                paste('High Expr  (N=',nH,')',sep='')),  
                      #legend.title='group',
                      xlab = 'Overall Survival (months)', ylab = 'Survival Probability',
                      #xlab = paste(type,'(months)'), ylab = 'Survival Probability',
                      font.x = c(11), font.y = c(11), ylim=c(0,1), #16
                      censor.size=1.5, size = 0.5,
                      ggtheme = theme_bw()+ theme(axis.line = element_line(colour = "black"),
                                                  #panel.grid.major = element_blank(),
                                                  #panel.grid.minor = element_blank(),
                                                  #panel.border = element_rect(colour='black'),
                                                  panel.border = element_blank(),
                                                  panel.background = element_blank(),
                                                  legend.text = element_text(size=12),#14
                                                  legend.title = element_blank(),
                                                  legend.position = 'none',
                                                  axis.text = element_text(size=9, color='black'))) #+
      
      ggplotly(p[[1]], height = 300, width = 325)
      
    })
    
    
    ##### highly expressed miRNAs
    output$high.expr.table.tcga <- DT::renderDataTable({
      
      expr.high.tcga[[project]]
      
      
    }, 
    callback=JS('$("button.buttons-copy").css("background","lightskyblue");
                 $("button.buttons-copy").css("border","white");
                 $("button.buttons-csv").css("background","lightskyblue");
                 $("button.buttons-csv").css("border","white");
                 $("button.buttons-excel").css("background","lightskyblue");
                 $("button.buttons-excel").css("border","white");
                 return table;'),
    options = list(pageLength = 10, dom = 'lfrtBp', buttons = c('copy', 'csv', 'excel')), #i
    extensions = "Buttons",
    selection = list(mode='none', selected=1), ### === not selectable
    rownames = FALSE,
    server = FALSE
    )
    
    
    output$high.expr.barplot.tcga <- renderPlot({
      
      dataForBarPlot <- expr.high.tcga[[project]][1:50,]
      
      dataForBarPlot$miRNA.ID <- factor(dataForBarPlot$miRNA.ID, levels=dataForBarPlot$miRNA.ID)
      
      p <- mirBarPlotFun(dataForBarPlot)
      
      #p <- ggplotly(p, tooltip=c("x", "y"))
      p
      
    })
    
    
    
    observeEvent(input$tcga_metadata, {
      
      meta.name <- input$tcga_metadata
      
      keep <- rowSums(expr > 0) >= 0.5*ncol(expr)
      expr.med <- apply(expr[keep,], 1, median)
      
      o <- order(expr.med, decreasing = T)
      expr.med <- expr.med[o]
      
      mir.id <- names(expr.med)
      mir.name <- mir.annotation[mir.id,]$Name
      
      
      if (meta.name=='sample_type') {
        expr.diy <- expr[mir.id,]
        groups.all <- meta[,meta.name]
        
      } else {
        samples <- which(meta$sample_type=='Tumor')
        expr.diy <- expr[mir.id,samples]
        
        if (meta.name=='age_at_initial_pathologic_diagnosis') {
          #cutoff <- reactiveVal(input$age)
          cutoff <- 65
          age <- meta[samples,meta.name]
          
          #groups.all <- ifelse(age>input$age, paste0('Age > ', input$age), paste0('Age <= ', input$age))
          groups.all <- ifelse(age>cutoff, paste0('Age > ', cutoff), paste0('Age <= ', cutoff))
          
        } else if (meta.name=='preop_psa' & project=='TCGA-PRAD') {
          #cutoff <- input$psa
          cutoff <- 20
          psa <- meta[samples,meta.name]
          
          groups.all <- ifelse(psa>cutoff, paste0('Preop PSA > ', cutoff), paste0('Preop PSA <= ', cutoff))
          
        }
        
        else {
          groups.all <- meta[samples,meta.name]
          
        }
      }
      
      filter <- which(rowSums(is.na(expr.diy))>0)
      if (length(filter)>0) {
        expr.diy <- expr.diy[-filter,]
      }
      
      groups.diy <- as.data.frame(table(groups.all), stringsAsFactors=F)
      
      group.names.diy <- sapply(groups.diy[,1], function(x) digest(x,algo='murmur32',seed=sample(1e9,1)))
      
      groups.diy <- cbind(group1=rep(NA, nrow(groups.diy)), group2=rep(NA, nrow(groups.diy)), groups.diy)
      
      groups.diy[,1] <- sprintf(
        '<input type="radio" name="%s" value="%s"/>',
        group.names.diy, 1)
      
      groups.diy[,2] <- sprintf(
        '<input type="radio" name="%s" value="%s"/>',
        group.names.diy, 2)
      
      
      colnames(groups.diy) <- c('Control','Case','Groups','N')
      
      output$groups.tcga.diy <- DT::renderDataTable(groups.diy, rownames = FALSE, escape = FALSE, selection = 'none', server = FALSE,
                                                    options=list(dom = 'tp', paging = TRUE, pageLength = 100, #ordering = FALSE,
                                                                 initComplete = JS("
                                                                                   function(setting, json) {
                                                                                   $(this.api().table().container())
                                                                                   .find('div.dataTables_paginate')
                                                                                   .css('display', this.api().page.info().pages <= 1 ? 'none' : 'block');
                                                                                   }"),
        drawCallback = JS("
                          function(settings) {
                          Shiny.unbindAll(this.api().table().node());
                          Shiny.bindAll(this.api().table().node());
                          }")
        ),
        callback = JS("
                      table.rows().every(function(i, tab, row) {
                      var $this = $(this.node());
                      //$(\"input[name='\" + this.data()[2] + \"']\").prop('checked', false);
                      //console.log($this.children()[0]);
                      //console.log($.parseHTML(this.data()[0])[0].name);
                      $this.attr('id', $.parseHTML(this.data()[0])[0].name); //one time hash value
                      //$this.attr('id', this.data()[2]); //Group Name
                      $this.addClass('shiny-input-radiogroup');
                      //console.log($this.prop('checked'));
                      });
                      Shiny.unbindAll(table.table().node());
                      Shiny.bindAll(table.table().node());
                      "
                      )
                      )
      
      
      shinyjs::hide('table_sample_type_tcga')
      shinyjs::hide('volcano_sample_type_tcga')
      
      observeEvent(input$deg.submit.tcga, {
        
        idx.diy <- unlist(sapply(group.names.diy, function(i) input[[i]]))
        
        req(length(unique(idx.diy))==2)
        
        groups <- names(idx.diy)
        
        idx1 <- which(idx.diy=='1')
        idx2 <- which(idx.diy=='2')
        
        control.groups <- groups[idx1]
        case.groups <- groups[idx2]
        
        deg.group <- groups.all
        
        idx <- which(deg.group %in% groups)
        deg.group <- deg.group[idx]
        
        
        deg.group.diy <- ifelse(deg.group %in% control.groups, 'Control', 'Case')
        deg.group.diy <- factor(deg.group.diy)
        
        design <- model.matrix(~0+deg.group.diy)
        colnames(design) <- levels(deg.group.diy)
        
        contrast.matrix <- makeContrasts(contrasts='Case - Control',
                                         levels=design)
        contrast.matrix
        
        ### Differential gene expression analysis (limma)
        
        fit <- lmFit(expr.diy[,idx], design)
        fit2 <- contrasts.fit(fit, contrast.matrix)
        fit2 <- eBayes(fit2)
        
        dgeTable <- topTable(fit2, coef=1, n=Inf, adjust.method='BH', sort.by='p')
        
        dataForVolcanoPlot <- dgeTable
        
        logFcThreshold <- log2(as.numeric(input$foldchange.tcga))
        adjPvalThreshold <- as.numeric(input$fdr.tcga)
        
        dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,
                                             logFC < logFcThreshold | adj.P.Val > adjPvalThreshold)] <- 'NS'
        dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,
                                             logFC >= logFcThreshold & adj.P.Val <= adjPvalThreshold)] <- 'UP'
        dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,
                                             logFC <= -logFcThreshold & adj.P.Val <= adjPvalThreshold)] <- 'DOWN'
        
        dataForVolcanoPlot <- data.frame(miRNA.Accession=rownames(dataForVolcanoPlot),
                                         miRNA.ID=mir.annotation[rownames(dataForVolcanoPlot),]$Name,
                                         dataForVolcanoPlot,
                                         stringsAsFactors = F)
        
        shinyjs::show('table_sample_type_tcga')
        shinyjs::show('volcano_sample_type_tcga')
        
        output$volcano_sample_type_tcga <- renderPlot({
          
          p <- volcanoPlotFun(dataForVolcanoPlot, logFcThreshold, adjPvalThreshold)
          p
          
        })
        
        
        output$table_sample_type_tcga <- DT::renderDataTable({
          
          dataForVolcanoPlot[,3:8] <- apply(dataForVolcanoPlot[,3:8], 2, 
                                            function(v) format(as.numeric(v), digits=3))
          dataForVolcanoPlot
          
        }, 
        callback=JS('$("button.buttons-copy").css("background","lightskyblue");
                 $("button.buttons-copy").css("border","white");
                 $("button.buttons-csv").css("background","lightskyblue");
                 $("button.buttons-csv").css("border","white");
                 $("button.buttons-excel").css("background","lightskyblue");
                 $("button.buttons-excel").css("border","white");
                 return table;'),
        options = list(pageLength = 10, dom = 'lfrtBp', buttons = c('copy', 'csv', 'excel')), #i
        extensions = "Buttons",
        selection = list(mode='none', selected=1), ### === not selectable
        rownames = FALSE,
        server = FALSE
        )
        
        
      })
      
    })
    
    ##### PCA analysis
    output$pca.tcga.2d <- renderPlotly({
      
      dataForPCAPlot <- pca.tcga[[project]]
      
      pc1 <- dataForPCAPlot$pc1[1]
      pc2 <- dataForPCAPlot$pc2[1]
      pc3 <- dataForPCAPlot$pc3[1]
      
      if (length(unique(dataForPCAPlot$Disease.Status))==2) {
        colors <- pie.colors[c(1,3)]
      } else if (unique(dataForPCAPlot$Disease.Status)=='Tumor') {
        colors <- pie.colors[3]
      } else if (unique(dataForPCAPlot$Disease.Status)=='Normal') {
        colors <- pie.colors[1]
      }
      
      
      p <- plot_ly(dataForPCAPlot, x = ~PC1, y = ~PC2,
                   color = ~Disease.Status,
                   colors = colors, 
                   alpha = 1,
                   marker = list(size=8),
                   text = ~Disease.Status,
                   hoverinfo = 'text',
                   showlegend=FALSE, height=400, width=500
                   
      )
      
      p <- p %>% layout(xaxis = list(title = paste0('PC1 (', pc1, '%)')),
                        yaxis = list(title = paste0('PC2 (', pc2, '%)')))
      
      p
      
    }) # }, height = 700, width = 700)
    
    output$pca.tcga.3d <- renderPlotly({
      
      dataForPCAPlot <- pca.tcga[[project]]
      
      pc1 <- dataForPCAPlot$pc1[1]
      pc2 <- dataForPCAPlot$pc2[1]
      pc3 <- dataForPCAPlot$pc3[1]
      
      if (length(unique(dataForPCAPlot$Disease.Status))==2) {
        colors <- pie.colors[c(1,3)]
      } else if (unique(dataForPCAPlot$Disease.Status)=='Tumor') {
        colors <- pie.colors[3]
      } else if (unique(dataForPCAPlot$Disease.Status)=='Normal') {
        colors <- pie.colors[1]
      }
      
      
      p <- plot_ly(dataForPCAPlot, x = ~PC1, y = ~PC2, z = ~PC3,
                   color = ~Disease.Status,
                   colors = colors, 
                   alpha = 1,
                   marker = list(size=5),
                   text = ~Disease.Status,
                   hoverinfo = 'text',
                   showlegend=FALSE, height=400, width=500
                   
      )
      
      p <- p %>% layout(scene = list(xaxis = list(title = paste0('PC1 (', pc1, '%)')),
                                     yaxis = list(title = paste0('PC2 (', pc2, '%)')),
                                     zaxis = list(title = paste0('PC3 (', pc3, '%)'))))
      
      p
      
      
    }) # }, height = 700, width = 700)
    
    
    ##### ROC analysis
    output$tcga_roc_analysis_table <- DT::renderDataTable({
      
      roc.tcga[[project]]
      
    },
    callback=JS('$("button.buttons-copy").css("background","lightskyblue");
                 $("button.buttons-copy").css("border","white");
                 $("button.buttons-csv").css("background","lightskyblue");
                 $("button.buttons-csv").css("border","white");
                 $("button.buttons-excel").css("background","lightskyblue");
                 $("button.buttons-excel").css("border","white");
                 return table;'),
    options = list(pageLength = 10, dom = 'lfrtBp', buttons = c('copy', 'csv', 'excel')), #i
    extensions = "Buttons",
    selection = list(mode='none', selected=1), ### === not selectable
    rownames = FALSE,
    server = FALSE
    )
    
    
    
    ##### feature selection
    
    output$tcga.feature.plot1 <- renderPlot({
      
      group <- meta$sample_type
      
      if (length(unique(group))==2 & (! project %in% c('TCGA-SKCM','TCGA-THYM','TCGA-CESC'))) {
        
        plot(tcga.feature.plot[[project]])
        
      } else {
        tcga.feature.plot[[project]]
      }
      
    }, height = 400)
    
    
    output$tcga.feature.table <- DT::renderDataTable({
      idx <- which(tcga.feature.table[[project]][,1]=='(Intercept)')
      if (length(idx)==1) {
        tcga.feature.table[[project]][-idx,]
      } else {
        tcga.feature.table[[project]]
      }
      
    },
    callback=JS('$("button.buttons-copy").css("background","lightskyblue");
                 $("button.buttons-copy").css("border","white");
                 $("button.buttons-csv").css("background","lightskyblue");
                 $("button.buttons-csv").css("border","white");
                 $("button.buttons-excel").css("background","lightskyblue");
                 $("button.buttons-excel").css("border","white");
                 return table;'),
    options = list(pageLength = 10, dom = 'lfrtBp', buttons = c('copy', 'csv', 'excel')), #i
    extensions = "Buttons",
    selection = list(mode='none', selected=1), ### === not selectable
    rownames = FALSE,
    server = FALSE
    )
    
    
    
    
    output$table_km_tcga <- DT::renderDataTable({
      
      km.tcga[[project]]
      
    },
    callback=JS('$("button.buttons-copy").css("background","lightskyblue");
                 $("button.buttons-copy").css("border","white");
                $("button.buttons-csv").css("background","lightskyblue");
                $("button.buttons-csv").css("border","white");
                $("button.buttons-excel").css("background","lightskyblue");
                $("button.buttons-excel").css("border","white");
                return table;'),
    options = list(pageLength = 10, dom = 'lfrtBp', buttons = c('copy', 'csv', 'excel')), #i
    extensions = "Buttons",
    selection = list(mode='none', selected=1), ### === not selectable
    rownames = FALSE,
    server = FALSE
    )
    
    output$table_coxph_tcga <- DT::renderDataTable({
      
      coxph.tcga[[project]]
      
    },
    callback=JS('$("button.buttons-copy").css("background","lightskyblue");
                 $("button.buttons-copy").css("border","white");
                $("button.buttons-csv").css("background","lightskyblue");
                $("button.buttons-csv").css("border","white");
                $("button.buttons-excel").css("background","lightskyblue");
                $("button.buttons-excel").css("border","white");
                return table;'),
    options = list(pageLength = 10, dom = 'lfrtBp', buttons = c('copy', 'csv', 'excel')), #i
    extensions = "Buttons",
    selection = list(mode='none', selected=1), ### === not selectable
    rownames = FALSE,
    server = FALSE
    )
    
    output$lasso_plot_tcga <- renderPlot({
      
      plot(lasso.plot.tcga[[project]])
      
      # if (project=='TCGA-GBM') {
      #   lasso.plot.tcga[[project]]
      # } else {
      #   plot(lasso.plot.tcga[[project]])
      # }
      
    })
    
    output$table_lasso_tcga <- DT::renderDataTable({
      
      lasso.tcga[[project]]
      
    },
    callback=JS('$("button.buttons-copy").css("background","lightskyblue");
                 $("button.buttons-copy").css("border","white");
                $("button.buttons-csv").css("background","lightskyblue");
                $("button.buttons-csv").css("border","white");
                $("button.buttons-excel").css("background","lightskyblue");
                $("button.buttons-excel").css("border","white");
                return table;'),
    options = list(pageLength = 10, dom = 'lfrtBp', buttons = c('copy', 'csv', 'excel')), #i
    extensions = "Buttons",
    selection = list(mode='none', selected=1), ### === not selectable
    rownames = FALSE,
    server = FALSE
    )
    
    output$risk_plot_tcga <- renderPlot({
      
      dataForKMPlot <- risk.km.data.tcga[[project]]
      
      if (is.null(dataForKMPlot)) {
        p <- ggplot(dataForKMPlot)
      } else {
        p <- KMRiskPlotFun(dataForKMPlot)
      }
      
      p
      
    })
    
    output$surv_roc_plot_tcga <- renderPlot({
      
      surv.roc.plot.tcga[[project]]
      
    })
    
    })
  
  
  
  #################### miRNomes
  
  output$ccma_datasets <- DT::renderDataTable({ccma.primary},
                                              options = list(pageLength = 5),
                                              selection = list(mode='single', selected=1)
  )
  
  observeEvent(input$ccma_datasets_rows_selected, {
    
    
    idx <- input$ccma_datasets_rows_selected
    req(idx)
    
    dataset <- as.character(ccma.datasets[idx,'Dataset'])
    
    meta <- meta.ccma[[dataset]]
    #expr <- expr.ccma[[dataset]]
    expr <- getCCMATable(dataset)
    
    groups <- meta$Group
    
    group.levels <- as.character(unlist(sapply(ccma.datasets[idx,'Group'], 
                                               function(x) strsplit(x, '; ')[[1]])))
    
    groups <- factor(groups, levels = group.levels)
    
    groups <- as.data.frame(table(groups), stringsAsFactors=F)
    
    disease.groups <- meta %>% group_by(Group) %>%
      summarise(Disease.Status=unique(Disease.Status))
    
    disease.status <- disease.groups$Disease.Status[match(groups$groups, disease.groups$Group)]
    
    ### deg
    deg.group.names <- sapply(groups[,1], function(x) digest(x,algo='murmur32',seed=sample(1e9,1)))
    deg.groups <- cbind(rep(NA, nrow(groups)), rep(NA, nrow(groups)), disease.status, groups)
    
    deg.groups[,1] <- sprintf(
      '<input type="radio" name="%s" value="%s"/>',
      deg.group.names, 1)
    
    deg.groups[,2] <- sprintf(
      '<input type="radio" name="%s" value="%s"/>',
      deg.group.names, 2)
    
    colnames(deg.groups) <- c('Control','Case','Disease.Status', 'Groups','N')
    
    
    
    ### roc
    
    roc.group.names <- sapply(groups[,1], function(x) digest(x,algo='murmur32',seed=sample(1e9,1)))
    roc.groups <- cbind(rep(NA, nrow(groups)), rep(NA, nrow(groups)), disease.status, groups)
    
    roc.groups[,1] <- sprintf(
      '<input type="radio" name="%s" value="%s"/>',
      roc.group.names, 1)
    
    roc.groups[,2] <- sprintf(
      '<input type="radio" name="%s" value="%s"/>',
      roc.group.names, 2)
    
    colnames(roc.groups) <- c('Control','Case','Disease.Status', 'Groups','N')
    
    
    ### feature selection
    
    feature.group.names <- sapply(groups[,1], function(x) digest(x,algo='murmur32',seed=sample(1e9,1)))
    feature.groups <- cbind(rep(NA, nrow(groups)), rep(NA, nrow(groups)), disease.status, groups)
    
    feature.groups[,1] <- sprintf(
      '<input type="radio" name="%s" value="%s"/>',
      feature.group.names, 1)
    
    feature.groups[,2] <- sprintf(
      '<input type="radio" name="%s" value="%s"/>',
      feature.group.names, 2)
    
    colnames(feature.groups) <- c('Control','Case','Disease.Status', 'Groups','N')
    
    ################
    output$groups.ccma <- DT::renderDataTable(deg.groups, rownames = FALSE, escape = FALSE, selection = 'none', server = FALSE,
                                              options=list(dom = 'tp', paging = TRUE, pageLength = 100, #ordering = FALSE,
                                                           initComplete = JS("
                                                                             function(setting, json) {
                                                                             $(this.api().table().container())
                                                                             .find('div.dataTables_paginate')
                                                                             .css('display', this.api().page.info().pages <= 1 ? 'none' : 'block');
                                                                             }"),
        drawCallback = JS("
                          function(settings) {
                          Shiny.unbindAll(this.api().table().node());
                          Shiny.bindAll(this.api().table().node());
                          }")
        ),
        callback = JS("
                      table.rows().every(function(i, tab, row) {
                      var $this = $(this.node());
                      //$(\"input[name='\" + this.data()[2] + \"']\").prop('checked', false);
                      //console.log($this.children()[0]);
                      //console.log($.parseHTML(this.data()[0])[0].name);
                      $this.attr('id', $.parseHTML(this.data()[0])[0].name); //one time hash value
                      //$this.attr('id', this.data()[2]); //Group Name
                      $this.addClass('shiny-input-radiogroup');
                      //console.log($this.prop('checked'));
                      });
                      Shiny.unbindAll(table.table().node());
                      Shiny.bindAll(table.table().node());
                      "
                      )
        
                      )
    
    output$roc.groups.ccma <- DT::renderDataTable(roc.groups, rownames = FALSE, escape = FALSE, selection = 'none', server = FALSE,
                                                  options=list(dom = 'tp', paging = TRUE, pageLength = 100, #ordering = FALSE,
                                                               initComplete = JS("
                                                                                 function(setting, json) {
                                                                                 $(this.api().table().container())
                                                                                 .find('div.dataTables_paginate')
                                                                                 .css('display', this.api().page.info().pages <= 1 ? 'none' : 'block');
                                                                                 }"),
        drawCallback = JS("
                          function(settings) {
                          Shiny.unbindAll(this.api().table().node());
                          Shiny.bindAll(this.api().table().node());
                          }")
        ),
        callback = JS("
                      table.rows().every(function(i, tab, row) {
                      var $this = $(this.node());
                      //$(\"input[name='\" + this.data()[2] + \"']\").prop('checked', false);
                      //console.log($this.children()[0]);
                      //console.log($.parseHTML(this.data()[0])[0].name);
                      $this.attr('id', $.parseHTML(this.data()[0])[0].name); //one time hash value
                      //$this.attr('id', this.data()[2]); //Group Name
                      $this.addClass('shiny-input-radiogroup');
                      //console.log($this.prop('checked'));
                      });
                      Shiny.unbindAll(table.table().node());
                      Shiny.bindAll(table.table().node());
                      "
                      )
        
                      )
    
    output$feature.selection.groups.ccma <- DT::renderDataTable(feature.groups, rownames = FALSE, escape = FALSE, selection = 'none', server = FALSE,
                                                                options=list(dom = 'tp', paging = TRUE, pageLength = 100, #ordering = FALSE,
                                                                             initComplete = JS("
                                                                                               function(setting, json) {
                                                                                               $(this.api().table().container())
                                                                                               .find('div.dataTables_paginate')
                                                                                               .css('display', this.api().page.info().pages <= 1 ? 'none' : 'block');
                                                                                               }"),
        drawCallback = JS("
                          function(settings) {
                          Shiny.unbindAll(this.api().table().node());
                          Shiny.bindAll(this.api().table().node());
                          }")
        ),
        callback = JS("
                      table.rows().every(function(i, tab, row) {
                      var $this = $(this.node());
                      //$(\"input[name='\" + this.data()[2] + \"']\").prop('checked', false);
                      //console.log($this.children()[0]);
                      //console.log($.parseHTML(this.data()[0])[0].name);
                      $this.attr('id', $.parseHTML(this.data()[0])[0].name); //one time hash value
                      //$this.attr('id', this.data()[2]); //Group Name
                      $this.addClass('shiny-input-radiogroup');
                      //console.log($this.prop('checked'));
                      });
                      Shiny.unbindAll(table.table().node());
                      Shiny.bindAll(table.table().node());
                      "
                      )
        
                      )
    
    
    
    
    output$dataset_summary <- renderText({ 
      dataset_summary <- as.character(paste0(ccma.datasets[idx,'Accession'], ': ', ccma.datasets[idx,'Title']))
      dataset_summary
    })
    
    
    output$gse <- renderUI({
      link <- as.character(ccma.datasets[idx,'Links'])
      tags$iframe(src=link, seamless="seamless", width='100%', height='600')
    })
    
    
    output$pie_disease_status <- renderPlotly({
      sample.freq <- table(meta$Disease.Status)
      dataForPiePlot <- data.frame(num=as.numeric(sample.freq), sam=names(sample.freq))
      
      p <- piePlotlyFun(dataForPiePlot)
      p
    })
    
    output$pie_group <- renderPlotly({
      sample.freq <- table(meta$Group)
      dataForPiePlot <- data.frame(num=as.numeric(sample.freq), sam=names(sample.freq))
      
      #p <- pieplotFun(dataForPiePlot)
      p <- piePlotlyFun(dataForPiePlot)
      p
    })
    
    
    output$high.expr.table.ccma <- DT::renderDataTable({
      
      expr.high.ccma[[dataset]]
      
    },
    callback=JS('$("button.buttons-copy").css("background","lightskyblue");
                 $("button.buttons-copy").css("border","white");
                $("button.buttons-csv").css("background","lightskyblue");
                $("button.buttons-csv").css("border","white");
                $("button.buttons-excel").css("background","lightskyblue");
                $("button.buttons-excel").css("border","white");
                return table;'),
    options = list(pageLength = 10, dom = 'lfrtBp', buttons = c('copy', 'csv', 'excel')), #i
    extensions = "Buttons",
    selection = list(mode='none', selected=1), ### === not selectable
    rownames = FALSE,
    server = FALSE
    )
    
    output$high.expr.barplot.ccma <- renderPlot({
      
      dataForBarPlot <- expr.high.ccma[[dataset]][1:50,]
      dataForBarPlot$miRNA.ID <- factor(dataForBarPlot$miRNA.ID, levels=dataForBarPlot$miRNA.ID)
      
      p <- mirBarPlotCCMAFun(dataForBarPlot)
      
      #p <- ggplotly(p, tooltip=c("x", "y"))
      p
      
    })
    
    shinyjs::hide('table_sample_type')
    shinyjs::hide('volcano_sample_type')
    
    observeEvent(input$deg.submit, {
      
      idx <- unlist(sapply(deg.group.names, function(i) input[[i]]))
      
      req(length(unique(idx))==2)
      
      groups <- names(idx)
      
      idx1 <- which(idx=='1')
      idx2 <- which(idx=='2')
      
      control.groups <- groups[idx1]
      case.groups <- groups[idx2]
      
      deg.group <- meta$Group
      
      idx <- which(deg.group %in% groups)
      deg.group <- deg.group[idx]
      
      
      deg.group.ccma <- ifelse(deg.group %in% control.groups, 'Control', 'Case')
      deg.group.ccma <- factor(deg.group.ccma)
      
      ###
      #req(length(levels(deg.group.ccma)==2))
      
      design <- model.matrix(~0+deg.group.ccma)
      colnames(design) <- levels(deg.group.ccma)
      
      contrast.matrix <- makeContrasts(contrasts='Case - Control',
                                       levels=design)
      contrast.matrix
      
      ### Differential gene expression analysis (limma)
      
      genes <- expr.high.ccma[[dataset]]$miRNA.Accession
      
      fit <- lmFit(expr[genes,idx], design)
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      
      dgeTable <- topTable(fit2, coef=1, n=Inf, adjust.method='BH', sort.by='p')
      
      dataForVolcanoPlot <- dgeTable
      
      logFcThreshold <- log2(input$foldchange)
      adjPvalThreshold <- input$fdr
      
      dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,
                                           logFC < logFcThreshold | adj.P.Val > adjPvalThreshold)] <- 'NS'
      dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,
                                           logFC >= logFcThreshold & adj.P.Val <= adjPvalThreshold)] <- 'UP'
      dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,
                                           logFC <= -logFcThreshold & adj.P.Val <= adjPvalThreshold)] <- 'DOWN'
      
      dataForVolcanoPlot <- data.frame(miRNA.Accession=rownames(dataForVolcanoPlot),
                                       miRNA.ID=mir.annotation[rownames(dataForVolcanoPlot),]$Name,
                                       dataForVolcanoPlot,
                                       stringsAsFactors = F)
      
      shinyjs::show('table_sample_type')
      shinyjs::show('volcano_sample_type')
      
      output$volcano_sample_type <- renderPlot({
        
        p <- volcanoPlotFun(dataForVolcanoPlot, logFcThreshold, adjPvalThreshold)
        p
        
        
      })
      
      
      output$table_sample_type <- DT::renderDataTable({
        
        dataForVolcanoPlot[,3:8] <- apply(dataForVolcanoPlot[,3:8], 2,
                                          function(v) format(as.numeric(v), digits=3))
        dataForVolcanoPlot
        
      }, 
      callback=JS('$("button.buttons-copy").css("background","lightskyblue");
                 $("button.buttons-copy").css("border","white");
                $("button.buttons-csv").css("background","lightskyblue");
                $("button.buttons-csv").css("border","white");
                $("button.buttons-excel").css("background","lightskyblue");
                $("button.buttons-excel").css("border","white");
                return table;'),
      options = list(pageLength = 10, dom = 'lfrtBp', buttons = c('copy', 'csv', 'excel')), #i
      extensions = "Buttons",
      selection = list(mode='none', selected=1), ### === not selectable
      rownames = FALSE,
      server = FALSE
      )
      
      
    })
    
    shinyjs::hide('ccma_roc_analysis_table')
    
    observeEvent(input$roc.submit, {
      
      idx <- unlist(sapply(roc.group.names, function(i) input[[i]]))
      
      req(length(unique(idx))==2)
      
      shinyjs::show('hide.roc')
      
      groups <- names(idx)
      
      idx1 <- which(idx=='1')
      idx2 <- which(idx=='2')
      
      control.groups <- groups[idx1]
      case.groups <- groups[idx2]
      
      deg.group <- meta$Group
      
      idx <- which(deg.group %in% groups)
      deg.group <- deg.group[idx]
      
      deg.group.ccma <- ifelse(deg.group %in% control.groups, 'Control', 'Case')
      deg.group.ccma <- factor(deg.group.ccma)
      group <- deg.group.ccma
      
      mir.id <- expr.high.ccma[[dataset]]$miRNA.Accession
      mir.name <- mir.annotation[mir.id,]$Name
      
      if (length(unique(group))!=2) {
        dataForForestPlot <- data.frame(matrix(ncol = 4, nrow = 0))
        colnames(dataForForestPlot) <- c('AUC','Lower95','Upper95','miRNA.ID') #,'P.Value',
      } else {
        
        dataForROCAnalysis <- expr[mir.id,idx]
        
        filter <- which(rowSums(is.na(dataForROCAnalysis))>0)
        
        if (length(filter)>0) {
          dataForROCAnalysis <- dataForROCAnalysis[-filter,]
        }
        
        dataForForestPlot <- apply(dataForROCAnalysis, 1, function(expr) ROCAnalysisFun(expr, group))
        dataForForestPlot <- t(apply(dataForForestPlot, 2, as.numeric))
        
        dataForForestPlot <- data.frame(miRNA.Accession=mir.id,
                                        miRNA.ID=mir.name,
                                        dataForForestPlot,
                                        row.names = NULL,
                                        stringsAsFactors = F)
        
        colnames(dataForForestPlot)[3:5] <- c('AUC','Lower95','Upper95') #,'P.Value'
        
        # o <- order(dataForForestPlot$P.Value, decreasing = F)
        # dataForForestPlot <- dataForForestPlot[o,]
        
        o <- order(dataForForestPlot$AUC, decreasing = T)
        dataForForestPlot <- dataForForestPlot[o,]
        
      }
      
      shinyjs::show('ccma_roc_analysis_table')
      
      ##### ROC analysis
      output$ccma_roc_analysis_table <- DT::renderDataTable({
        
        dataForForestPlot
        
      },
      callback=JS('$("button.buttons-copy").css("background","lightskyblue");
                 $("button.buttons-copy").css("border","white");
                  $("button.buttons-csv").css("background","lightskyblue");
                  $("button.buttons-csv").css("border","white");
                  $("button.buttons-excel").css("background","lightskyblue");
                  $("button.buttons-excel").css("border","white");
                  return table;'),
      options = list(pageLength = 10, dom = 'lfrtBp', buttons = c('copy', 'csv', 'excel')), #i
      extensions = "Buttons",
      selection = list(mode='none', selected=1), ### === not selectable
      rownames = FALSE,
      server = FALSE
      )
      
    })
    
    shinyjs::hide('ccma.feature.plot1')
    shinyjs::hide('ccma.feature.table')
    
    observeEvent(input$feature.selection.submit, {
      
      idx <- unlist(sapply(feature.group.names, function(i) input[[i]]))
      
      req(length(unique(idx))==2)
      
      shinyjs::show('hide.feature')
      
      groups <- names(idx)
      
      idx1 <- which(idx=='1')
      idx2 <- which(idx=='2')
      
      control.groups <- groups[idx1]
      case.groups <- groups[idx2]
      
      deg.group <- meta$Group
      
      idx <- which(deg.group %in% groups)
      deg.group <- deg.group[idx]
      
      deg.group.ccma <- ifelse(deg.group %in% control.groups, 'Control', 'Case')
      deg.group.ccma <- factor(deg.group.ccma)
      group <- deg.group.ccma
      
      mir.id <- expr.high.ccma[[dataset]]$miRNA.Accession
      mir.name <- mir.annotation[mir.id,]$Name
      
      if (length(unique(group))!=2) {
        dataForForestPlot <- data.frame(matrix(ncol = 5, nrow = 0))
        colnames(dataForForestPlot) <- c('AUC','Lower95','Upper95','P.Value','miRNA.ID')
      } else {
        
        dataForROCAnalysis <- expr[mir.id,idx]
        
        filter <- which(rowSums(is.na(dataForROCAnalysis))>0)
        
        if (length(filter)>0) {
          dataForROCAnalysis <- dataForROCAnalysis[-filter,]
        }
        
        set.seed(777)
        cvfit<-cv.glmnet(x=t(dataForROCAnalysis),y=group, alpha=1, 
                         family = "binomial", type.measure="class") #class
        
      }
      
      shinyjs::show('ccma.feature.plot1')
      shinyjs::show('ccma.feature.table')
      
      ##### feature selection
      
      output$ccma.feature.plot1 <- renderPlot({
        
        plot(cvfit)
        
      }, height = 400)
      
      
      output$ccma.feature.table <- DT::renderDataTable({
        
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
        
        feature <- feature[-1,]
        
        o <- order(abs(feature$Coefficient), decreasing = T)
        feature <- feature[o,]
        
        feature[,3] <- round(feature[,3],4)
        
        feature
        
      },
      callback=JS('$("button.buttons-copy").css("background","lightskyblue");
                 $("button.buttons-copy").css("border","white");
                  $("button.buttons-csv").css("background","lightskyblue");
                  $("button.buttons-csv").css("border","white");
                  $("button.buttons-excel").css("background","lightskyblue");
                  $("button.buttons-excel").css("border","white");
                  return table;'),
      options = list(pageLength = 10, dom = 'lfrtBp', buttons = c('copy', 'csv', 'excel')), #i
      extensions = "Buttons",
      selection = list(mode='none', selected=1), ### === not selectable
      rownames = FALSE,
      server = FALSE
      )
      
    })
    
    
    output$pca.ccma.2d <- renderPlotly({
      
      dataForPCAPlot <- pca.ccma[[dataset]]
      
      pc1 <- dataForPCAPlot$pc1[1]
      pc2 <- dataForPCAPlot$pc2[1]
      pc3 <- dataForPCAPlot$pc3[1]
      
      colors <- pie.colors[-2][1:length(unique(dataForPCAPlot$Disease.Status))]
      
      p <- plot_ly(dataForPCAPlot, x = ~PC1, y = ~PC2,
                   color = ~Disease.Status,
                   colors = colors, 
                   alpha = 1,
                   marker = list(size=8),
                   text = ~Disease.Status,
                   hoverinfo = 'text',
                   showlegend=FALSE, height=400, width=500
                   
      )
      
      p <- p %>% layout(xaxis = list(title = paste0('PC1 (', pc1, '%)')),
                        yaxis = list(title = paste0('PC2 (', pc2, '%)')))
      
      p
      
      
    }) # }, height = 700, width = 700)
    
    
    output$pca.ccma.3d <- renderPlotly({
      
      dataForPCAPlot <- pca.ccma[[dataset]]
      
      pc1 <- dataForPCAPlot$pc1[1]
      pc2 <- dataForPCAPlot$pc2[1]
      pc3 <- dataForPCAPlot$pc3[1]
      
      colors <- pie.colors[-2][1:length(unique(dataForPCAPlot$Disease.Status))]
      
      
      p <- plot_ly(dataForPCAPlot, x = ~PC1, y = ~PC2, z = ~PC3,
                   color = ~Disease.Status,
                   colors = colors, 
                   alpha = 1,
                   marker = list(size=5),
                   text = ~Disease.Status,
                   hoverinfo = 'text',
                   showlegend=FALSE, height=400, width=500
                   
      )
      
      p <- p %>% layout(scene = list(xaxis = list(title = paste0('PC1 (', pc1, '%)')),
                                     yaxis = list(title = paste0('PC2 (', pc2, '%)')),
                                     zaxis = list(title = paste0('PC3 (', pc3, '%)'))))
      
      
      p
      
      
    }) # }, height = 700, width = 700)
    
    
    output$panelStatus <- reactive({
      idx %in% c(4,16:18,21:26,28,31:33,26,39:40)
    })
    outputOptions(output, "panelStatus", suspendWhenHidden = FALSE)
    
    
    
    output$pca.ccma.2d.group <- renderPlotly({
      
      dataForPCAPlot <- pca.ccma[[dataset]]
      
      pc1 <- dataForPCAPlot$pc1[1]
      pc2 <- dataForPCAPlot$pc2[1]
      pc3 <- dataForPCAPlot$pc3[1]
      
      colors <- pie.colors[-2][1:length(unique(dataForPCAPlot$Group))]
      
      p <- plot_ly(dataForPCAPlot, x = ~PC1, y = ~PC2,
                   color = ~Group,
                   colors = colors, 
                   alpha = 1,
                   marker = list(size=8),
                   text = ~Group,
                   hoverinfo = 'text',
                   showlegend=FALSE, height=400, width=500
                   
      )
      
      p <- p %>% layout(xaxis = list(title = paste0('PC1 (', pc1, '%)')),
                        yaxis = list(title = paste0('PC2 (', pc2, '%)')))
      
      p
      
      
    }) # }, height = 700, width = 700)
    
    output$pca.ccma.3d.group <- renderPlotly({
      
      dataForPCAPlot <- pca.ccma[[dataset]]
      
      pc1 <- dataForPCAPlot$pc1[1]
      pc2 <- dataForPCAPlot$pc2[1]
      pc3 <- dataForPCAPlot$pc3[1]
      
      colors <- pie.colors[-2][1:length(unique(dataForPCAPlot$Group))]
      
      p <- plot_ly(dataForPCAPlot, x = ~PC1, y = ~PC2, z = ~PC3,
                   color = ~Group,
                   colors = colors, 
                   alpha = 1,
                   marker = list(size=5),
                   text = ~Group,
                   hoverinfo = 'text',
                   showlegend=FALSE, height=400, width=500
                   
      )
      
      p <- p %>% layout(scene = list(xaxis = list(title = paste0('PC1 (', pc1, '%)')),
                                     yaxis = list(title = paste0('PC2 (', pc2, '%)')),
                                     zaxis = list(title = paste0('PC3 (', pc3, '%)'))))
      
      
      
      p
      
      
    }) # }, height = 700, width = 700)
    
    
    
  })
  
  }